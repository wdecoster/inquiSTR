//! Thread-local BAM reader pool for efficient parallel processing
//!
//! This module provides a pool of BAM readers that can be used with rayon's parallel
//! iterators. Each thread gets its own BAM reader, avoiding the need for a global lock
//! while respecting file descriptor limits.

use log::debug;
use rust_htslib::bam::IndexedReader;
use std::cell::RefCell;
use std::sync::Mutex;

use crate::bam_utils::{get_bam_reader, get_bam_reader_with_validation};
use crate::errors::InquiSTRResult;

/// A pool of BAM readers, one per thread
///
/// Uses thread-local storage so each rayon thread gets its own reader.
/// Readers are lazily initialized on first use and reused for subsequent batches.
pub struct BamReaderPool {
    bam_path: String,
    reference: Option<String>,
    validate_phasing: bool,
    /// Mutex to serialize BAM reader creation (htslib safety requirement)
    creation_lock: Mutex<()>,
}

impl BamReaderPool {
    /// Create a new pool with the given BAM path and optional reference
    ///
    /// The first reader created will validate phasing if `validate_phasing` is true.
    pub fn new(bam_path: String, reference: Option<String>, validate_phasing: bool) -> Self {
        Self { bam_path, reference, validate_phasing, creation_lock: Mutex::new(()) }
    }

    /// Get or create a BAM reader for the current thread
    ///
    /// This uses thread-local storage via rayon's current thread index.
    /// The reader is created lazily on first access and reused for all subsequent calls
    /// from the same thread.
    fn create_reader(&self) -> InquiSTRResult<IndexedReader> {
        // Serialize reader creation to prevent concurrent htslib index operations
        let _guard = self.creation_lock.lock().unwrap();

        debug!("Creating new BAM reader for thread");

        if self.validate_phasing {
            get_bam_reader_with_validation(&self.bam_path, &self.reference)
        } else {
            get_bam_reader(&self.bam_path, &self.reference)
        }
    }

    /// Execute a closure with a BAM reader for the current thread
    ///
    /// The reader is obtained from thread-local storage (or created if first use).
    /// This is the main entry point for using the pool in parallel code.
    pub fn with_reader<F, T>(&self, f: F) -> InquiSTRResult<T>
    where
        F: FnOnce(&mut IndexedReader) -> InquiSTRResult<T>,
    {
        // Use thread-local storage for the reader
        // Each rayon thread will have its own reader
        thread_local! {
            static THREAD_READER: RefCell<Option<IndexedReader>> = const { RefCell::new(None) };
        }

        THREAD_READER.with(|cell| {
            let mut reader_opt = cell.borrow_mut();

            // Create reader if this thread doesn't have one yet
            if reader_opt.is_none() {
                *reader_opt = Some(self.create_reader()?);
            }

            // Use the reader
            f(reader_opt.as_mut().unwrap())
        })
    }
}

// Note: BamReaderPool is Sync because:
// - bam_path and reference are immutable String/Option<String>
// - validate_phasing is an immutable bool
// - creation_lock is Mutex<()> which is Sync
// Thread-local storage handles the actual readers safely
unsafe impl Sync for BamReaderPool {}

/// Pre-create BAM readers for a fixed number of threads
///
/// This is useful when you want to ensure all readers are opened before
/// parallel processing starts, to catch errors early and avoid opening
/// readers during processing.
pub fn create_readers_for_threads(
    bam_path: &str,
    reference: &Option<String>,
    num_threads: usize,
    validate_first: bool,
) -> InquiSTRResult<Vec<IndexedReader>> {
    // Serialize all reader creation
    let mut readers = Vec::with_capacity(num_threads);

    for i in 0..num_threads {
        debug!("Pre-creating BAM reader {} of {}", i + 1, num_threads);
        let reader = if i == 0 && validate_first {
            get_bam_reader_with_validation(&bam_path.to_string(), reference)?
        } else {
            get_bam_reader(&bam_path.to_string(), reference)?
        };
        readers.push(reader);
    }

    Ok(readers)
}
