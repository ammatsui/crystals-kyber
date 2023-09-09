use crate::{params::*, sample::*, poly::*, packing::*};
use rand::Rng;

/* CPAPKE procedures */

pub fn keyGen() -> ([u8; KEMPK_BYTES], [u8; KEMSK_BYTES])
{
    let mut sk = [0u8; KEMSK_BYTES];

    let (pk, skp) = crate::cpapke::keyGen();

    let z = rand::thread_rng().gen::<[u8; SEEDBYTES]>();
    
    let mut hpk = [0u8; SYMBYTES];
    H(&pk, &mut hpk);

    sk[..SK_BYTES].copy_from_slice(&skp[..SK_BYTES]);
    sk[SK_BYTES..PK_BYTES+SK_BYTES].copy_from_slice(&pk[..PK_BYTES]);
    sk[PK_BYTES+SK_BYTES..PK_BYTES+SK_BYTES+SYMBYTES].copy_from_slice(&hpk[..SYMBYTES]);
    sk[PK_BYTES+SK_BYTES+SYMBYTES..].copy_from_slice(&z[..SEEDBYTES]);
           
    (pk, sk)
}


pub fn encapsulation(pk: &[u8]) -> ([u8; CIPHERTEXTBYTES], [u8; SHAREDKEY])
{
    let seed = rand::thread_rng().gen::<[u8; SEEDBYTES]>();
    let mut m = [0u8; SEEDBYTES]; 
    H(&seed, &mut m);

    let mut tmp = [0u8; SYMBYTES];
    H(&pk, &mut tmp); 

    let mut seed = [0u8; SEEDBYTES+SYMBYTES]; 
    seed[..SEEDBYTES].copy_from_slice(&m);
    seed[SEEDBYTES..].copy_from_slice(&tmp);
    
    let mut res = [0u8; 2*SEEDBYTES];
    G(&seed, &mut res);
    let (mut shk, r) = res.split_at(SEEDBYTES); 

    let c = crate::cpapke::encryption(&pk, &m, &r);

    let mut tmp = [0u8; SYMBYTES];
    H(&c, &mut tmp); 

    let mut seed = [0u8; 2*SEEDBYTES]; 
    seed[..SEEDBYTES].copy_from_slice(&shk);
    seed[SEEDBYTES..].copy_from_slice(&tmp);

    let mut shk = [0u8; SHAREDKEY];
    kdf(&seed, &mut shk);

    (c, shk)

}


pub fn decapsulation(sk: &[u8], c: &[u8]) -> [u8; SHAREDKEY]
{
    let mut pk = [0u8; KEMPK_BYTES];
    pk[..KEMPK_BYTES].copy_from_slice(&sk[SK_BYTES..SK_BYTES+KEMPK_BYTES]);

    let mut h = [0u8; SYMBYTES];
    h[..SYMBYTES].copy_from_slice(&sk[SK_BYTES+KEMPK_BYTES..SK_BYTES+KEMPK_BYTES+SYMBYTES]);

    let mut z = [0u8; SEEDBYTES];
    z[..SEEDBYTES].copy_from_slice(&sk[SK_BYTES+KEMPK_BYTES+SYMBYTES..]);

    let mp = crate::cpapke::decryption(&sk, &c);

    let mut seed = [0u8; SEEDBYTES+SYMBYTES]; 
    seed[..SEEDBYTES].copy_from_slice(&mp);
    seed[SEEDBYTES..].copy_from_slice(&h);
    
    let mut res = [0u8; 2*SYMBYTES];
    G(&seed, &mut res);
    let (mut shkp, r) = res.split_at(SEEDBYTES);

    let cp = crate::cpapke::encryption(&pk, &mp, &r);

    let mut hc = [0u8; SYMBYTES];
    H(&c, &mut hc);

    let mut shk = [0u8; SEEDBYTES];
    if c == cp 
    {
        let mut seed = [0u8; SEEDBYTES+SYMBYTES]; 
        seed[..SEEDBYTES].copy_from_slice(&shkp);
        seed[SEEDBYTES..].copy_from_slice(&hc);
        kdf(&seed, &mut shk);
    }
    else  
    { 
        let mut seed = [0u8; SEEDBYTES+SYMBYTES]; 
        seed[..SEEDBYTES].copy_from_slice(&z);
        seed[SEEDBYTES..].copy_from_slice(&hc);
        kdf(&seed, &mut shk);
    }
    shk

}