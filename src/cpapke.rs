use crate::{params::*, sample::*, poly::*, packing::*};
use rand::Rng;

/* CPAPKE procedures */

pub fn keyGen() -> ([u8; PK_BYTES], [u8; SK_BYTES])
{
    let mut sk = [0u8; SK_BYTES];
    let mut pk = [0u8; PK_BYTES];

    let d  = rand::thread_rng().gen::<[u8; SEEDBYTES]>();
    

    let mut tmp = [0u8; 2*SEEDBYTES];
    G(&d, &mut tmp);
    let (rho, sigma) = tmp.split_at(SEEDBYTES);
   

    let A = get_matrix(&rho);
    let mut s = get_noise(&sigma);
    let mut e = get_error(&sigma);

    s = Ntt(&s);
    e = Ntt(&e);

    let mut t: VecPoly<K> = m_mult_v(&A, &s);
    t = v_add(&t, &e);

    t = modq(&t);
    s = modq(&s);

    pack_pk(&mut pk, &t, &rho);
    pack_sk(&mut sk, &s);
    
    (pk, sk)
}


pub fn encryption(pk: &[u8], m: &[u8], rc: &[u8]) -> [u8; CIPHERTEXTBYTES]
{
    let mut c = [0u8; CIPHERTEXTBYTES];
    let mut t = VecPoly::<K>::default();
    let mut rho = [0u8; SEEDBYTES];

    unpack_pk(&pk, &mut t, &mut rho);

    let mut A = get_matrix(&mut rho);
    

    let mut n = 0;
    let mut r = VecPoly::<K>::default();
    for i in 0..K
    {
        let mut tmp = [0u8; 2*SEEDBYTES*ETA1];
        let mut seed = [0u8; SEEDBYTES+2];
        seed[..SEEDBYTES].copy_from_slice(&rc[..SEEDBYTES]);
        seed[SEEDBYTES..].copy_from_slice(&[n as u8, (n>>8) as u8]);
        prf(&seed, &mut tmp);
        r.poly[i] = CBD(&tmp, ETA1);
        n+=1;
    }
    
    let mut e1 = VecPoly::<K>::default();
    for i in 0..K
    {
        let mut tmp = [0u8; 2*SEEDBYTES*ETA2];
        let mut seed = [0u8; SEEDBYTES+2];
        seed[..SEEDBYTES].copy_from_slice(&rc[..SEEDBYTES]);
        seed[SEEDBYTES..].copy_from_slice(&[n as u8, (n>>8) as u8]);
        prf(&seed, &mut tmp);
        e1.poly[i] = CBD(&tmp, ETA2);
        n+=1;
    }
    
    let mut tmp = [0u8; 2*SEEDBYTES*ETA2];
    let mut seed = [0u8; SEEDBYTES+2];
    seed[..SEEDBYTES].copy_from_slice(&rc[..SEEDBYTES]);
    seed[SEEDBYTES..].copy_from_slice(&[n as u8, (n>>8) as u8]);
    prf(&seed, &mut tmp);
    let e2 = CBD(&tmp, ETA2);
    
    /* Atr = rA */
    r = Ntt(&r);
    let mut u = inv_Ntt(&v_mult_m(&r, &A));
    u = v_add(&u, &e1);
    u = modq(&u);

    let mut v = inv_ntt(&v_mult_v(&r, &t)); 
    add(&mut v, &e2);

    let mut tmp = Poly::default();
    decode(&mut tmp, m, 1);
    tmp = decompress(&tmp, 1);
    
    add(&mut v, &tmp);
    v = _modq(&v);
    
    
    let mut tmp = Compress(&u, Du as i16);
    Encode(&mut c, &tmp, Du);
    let mut tmp = compress(&v, Dv as i16);
    encode(&mut c[Du*K*N/8..], &tmp, Dv);

    c
}


pub fn decryption(sk: &[u8], c: &[u8]) -> [u8; MESSAGEBYTES]
{
    let mut m = [0u8; MESSAGEBYTES];

    let mut tmp = VecPoly::<K>::default();
    Decode(&mut tmp, &c[..Du*K*N/8], Du);
    let u = Decompress(&mut tmp, Du as i16);

    let mut tmp = Poly::default();
    decode(&mut tmp, &c[Du*K*N/8..], Dv);
    let mut v = decompress(&mut tmp, Dv as i16);

    let mut s = VecPoly::<K>::default();
    unpack_sk(sk, &mut s);

    /* sTu = us */
    let us = inv_ntt(&v_mult_v(&Ntt(&u),&s));
    let us = neg(&us);
    add(&mut v, &us);

    v = _modq(&v);
    let tmp = compress(&v, 1);
    encode(&mut m, &tmp, 1);

    m
}