use rand::{self, RngExt, SeedableRng};
use std::thread;
use std::sync::{Mutex, Arc, Condvar};

type WorkItem<R> = Box<dyn FnOnce(&mut R) -> () + Send>;

pub struct ThreadPool<R : RngExt + SeedableRng + Send> {
  work_queue : Arc<Mutex<Vec<WorkItem<R>>>>,
  condvar : Arc<Condvar>,
}

impl<R : RngExt + SeedableRng + Send + 'static> ThreadPool<R> {
  pub fn new(rng : &mut R, thread_count : usize) -> Self {
    let work_queue : Arc<Mutex<Vec<WorkItem<R>>>> = Arc::new(Mutex::new(vec![]));
    let condvar = Arc::new(Condvar::new());
    for _ in 0..thread_count {
      let thread_work_queue = Arc::clone(&work_queue);
      let thread_condvar = Arc::clone(&condvar);
      let mut thread_rng = rng.fork();
      thread::spawn(move || {
        loop {
          let mut mutex_guard = thread_work_queue.lock().unwrap();
          if let Some(work) = mutex_guard.pop() {
            drop(mutex_guard);
            work(&mut thread_rng);
          } else {
            let new_lock = thread_condvar.wait(mutex_guard).unwrap();
            drop(new_lock);
          }
        }
      });
    }
    ThreadPool {work_queue, condvar}
  }
  
  pub fn submit<F : FnOnce(&mut R) -> () + Send + 'static>(&self, work : F) {
    self.work_queue.lock().unwrap().push(Box::new(work));
    self.condvar.notify_all();
  }
}
