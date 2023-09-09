use crate::{params::*, sample::*, poly::*, packing::*};
use rand::Rng;

/* CPAPKE procedures */

pub fn keyGen() -> ([u8; KEMPK_BYTES], [u8; KEMSK_BYTES])
{
    let mut sk = [0u8; KEMSK_BYTES];

    let (pk, skp) = crate::cpapke::keyGen();

    let z = rand::thread_rng().gen::<[u8; 32]>();
    
    let mut tmp = [0u8; 64];
    H(&pk, &mut tmp);

    sk[..SK_BYTES].copy_from_slice(&skp[..SK_BYTES]);
    sk[SK_BYTES..2*SK_BYTES].copy_from_slice(&pk[..SK_BYTES]);
    sk[2*SK_BYTES..2*SK_BYTES+64].copy_from_slice(&tmp[..64]);
    sk[2*SK_BYTES+64..].copy_from_slice(&z[..32]);
           
    (pk, sk)
}


pub fn encapsulation(pk: &[u8]) -> ([u8; CIPHERTEXTBYTES], [u8; SHAREDKEY])
{
    let seed = rand::thread_rng().gen::<[u8; 32]>();
    let mut m = [0u8; 64]; //size?
    H(&seed, &mut m);

    let mut tmp = [0u8; 64];
    H(&pk, &mut tmp); //size

    let mut seed = [0u8; 64+64]; //m+h(pk)
    seed[..64].copy_from_slice(&m);
    seed[64..].copy_from_slice(&tmp);
    
    let mut res = [0u8; 64];
    G(&seed, &mut res);
    let (mut shk, r) = res.split_at(32); //shk = shared key

    let c = crate::cpapke::encryption(&pk, &m, &r);

    let mut tmp = [0u8; 64];
    H(&c, &mut tmp); //size

    let mut seed = [0u8; 64]; //shk+h(c)
    seed[..32].copy_from_slice(&shk);
    seed[32..].copy_from_slice(&tmp);

    let mut shk = [0u8; SHAREDKEY];
    kdf(&seed, &mut shk);

    (c, shk)

}


pub fn decapsulation(sk: &[u8], c: &[u8]) -> [u8; SHAREDKEY]
{
    let mut pk = [0u8; KEMPK_BYTES];
    pk[..KEMPK_BYTES].copy_from_slice(&sk[SK_BYTES..PK_BYTES]);

    let mut h = [0u8; 64];
    h[..64].copy_from_slice(&sk[PK_BYTES+SK_BYTES..64]);

    let mut z = [0u8; 32];
    z[..32].copy_from_slice(&sk[PK_BYTES+SK_BYTES+64..]);

    let mp = crate::cpapke::decryption(&sk, &c);

    let mut seed = [0u8; 64+64]; //m+h
    seed[..64].copy_from_slice(&mp);
    seed[64..].copy_from_slice(&h);
    
    let mut res = [0u8; 64];
    G(&seed, &mut res);
    let (mut shkp, r) = res.split_at(32);

    let cp = crate::cpapke::encryption(&pk, &mp, &r);

    let mut hc = [0u8; 64];
    H(&c, &mut hc);

    let mut shk = [0u8; 32];
    if c == cp 
    {
        let mut seed = [0u8; 32+64]; //m+h
        seed[..32].copy_from_slice(&shkp);
        seed[32..].copy_from_slice(&hc);
        kdf(&seed, &mut shk);
    }
    else  
    { 
        let mut seed = [0u8; 32+64]; //m+h
        seed[..32].copy_from_slice(&z);
        seed[32..].copy_from_slice(&hc);
        kdf(&seed, &mut shk);
    }
    shk

}