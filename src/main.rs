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
    let m = b"hbvpedfghgjkhljjnghhhhhhonriobej";

//    let mut a = poly::Poly::default();
//    let mut aa = poly::Poly::default();
//    packing::decode(&mut a, m, 1);

//    let mut tmp = [0u8; 32];//SK_BYTES];
//    packing::encode(&mut tmp, a, 1);
//    assert_eq!(*m, tmp[..32]);
//    packing::decode(&mut aa, &tmp, 1);
//    assert_eq!(a.coeff, aa.coeff);

    //check

    let (pk, sk) = kex::keyGen();


    let rc  = rand::thread_rng().gen::<[u8; 32]>();
    let c = kex::encryption(&pk, m, &rc, &sk);

    let mm = kex::decryption(&sk, &c);

    assert_eq!(*m, mm);

  
}
