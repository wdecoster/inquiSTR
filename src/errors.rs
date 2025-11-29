/// Custom error types for inquiSTR
use std::fmt;

#[derive(Debug, Clone)]
pub struct InquiSTRError {
    pub message: String,
}

impl InquiSTRError {
    pub fn new(message: String) -> Self {
        InquiSTRError { message }
    }
}

impl fmt::Display for InquiSTRError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for InquiSTRError {}

impl From<std::io::Error> for InquiSTRError {
    fn from(err: std::io::Error) -> Self {
        InquiSTRError::new(format!("IO error: {}", err))
    }
}

impl From<String> for InquiSTRError {
    fn from(err: String) -> Self {
        InquiSTRError::new(err)
    }
}

impl From<&str> for InquiSTRError {
    fn from(err: &str) -> Self {
        InquiSTRError::new(err.to_string())
    }
}

pub type InquiSTRResult<T> = Result<T, InquiSTRError>;
