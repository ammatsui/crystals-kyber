pub mod params;
pub mod ntt;
pub mod utils;
pub mod poly;
pub mod sample;
pub mod packing;
pub mod kex;
pub mod kem;

use rand::Rng;

fn main()
{
    let (pk, sk) = kex::keyGen();

    let m = b"hbvpedfghgjkhljjnghhhhhhonriobej";

    let rc  = rand::thread_rng().gen::<[u8; 32]>();
    let c = kex::encryption(&pk, m, &rc);

    let mm = kex::decryption(&sk, &c);

    assert_eq!(*m, mm);

  
}
