#[cfg(feature = "mode512")]
mod mode512;
#[cfg(not(any(feature = "mode512", feature = "mode1024")))]
mod mode768;
#[cfg(feature = "mode1024")]
mod mode1024;

#[cfg(feature = "mode512")]
pub use mode512::*;
#[cfg(not(any(feature = "mode512", feature = "mode1024")))]
pub use mode768::*;
#[cfg(feature = "mode1024")]
pub use mode1024::*;


pub const Q: usize = 3329; 
pub const N: usize = 256;

pub const MESSAGEBYTES: usize = 32; 
pub const CIPHERTEXTBYTES: usize = Du*K*N/8 + Dv*N/8;


pub const SYMBYTES: usize = 32; //64
pub const SEEDBYTES: usize = 32;
pub const PACKED_KEYS: usize = 12; 

/* shared key */
pub const SHAREDKEY: usize = 32;

/* Kyber public key */
pub const PK_BYTES: usize = PACKED_KEYS*K*N/8 + SEEDBYTES;
/* Kyber secret key */
pub const SK_BYTES: usize = PACKED_KEYS*K*N/8;

/* encapsulation sizes */
pub const KEMSK_BYTES: usize = PK_BYTES+SK_BYTES+SEEDBYTES+SYMBYTES;//2*PACKED_KEYS*K*N/8 + SEEDBYTES + SYMBYTES;
pub const KEMPK_BYTES: usize = PK_BYTES;//PACKED_KEYS*K*N/8 + SEEDBYTES;
