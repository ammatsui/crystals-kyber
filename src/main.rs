pub mod params;
pub mod ntt;
pub mod utils;
pub mod poly;
pub mod sample;
pub mod packing;
pub mod cpapke;
pub mod ccakem;

use rand::Rng;

fn main()
{
    let m = b"hbvpedfghgjkhljjnghhhhhhonriobej";


    let (pk, sk) = cpapke::keyGen();


    let rc  = rand::thread_rng().gen::<[u8; 32]>();
    let c = cpapke::encryption(&pk, m, &rc);

    let mm = cpapke::decryption(&sk, &c);

    assert_eq!(*m, mm);

    let (pk, sk) = ccakem::keyGen();

    let (c, shk) = ccakem::encapsulation(&pk);

    let shkp = ccakem::decapsulation(&sk, &c);

    assert_eq!(shkp, shk);
}
