use crate::{params::*, sample::*, poly::*, packing::*};
use rand::Rng;

/* CPAPKE procedures */

pub fn keyGen() -> ([u8; PK_BYTES], [u8; SK_BYTES])
{
    let mut sk = [0u8; SK_BYTES];
    let mut pk = [0u8; PK_BYTES];

    let d  = rand::thread_rng().gen::<[u8; 32]>();
    
    let (mut rho, mut sigma) = ([0u8; 32], [0u8; 32]);

    let mut tmp = [0u8; 64];
    G(&d, &mut tmp);
    rho.copy_from_slice(&tmp[..32]);
    sigma.copy_from_slice(&tmp[32..]);

    let A = get_matrix(&rho);
    let mut s = get_noise(&sigma);
    let mut e = get_error(&sigma);

    s = Ntt(&s);
    e = Ntt(&e);

    let mut t: VecPoly<K> = m_mult_v(&A, &Ntt(&s));
    t = v_add(&t, &e);

    pack_pk(&mut pk, &t, &rho);
    pack_sk(&mut sk, &s);

    (pk, sk)
}


pub fn encryption(pk: &[u8], m: &[u8], rc: &[u8]) -> [u8; CIPHERTEXTBYTES]
{
    let mut c = [0u8; CIPHERTEXTBYTES];
    let mut t = VecPoly::<K>::default();
    let mut rho = [0u8; 32];

    unpack_pk(&pk, &mut t, &mut rho);

    let mut A = get_matrix(&mut rho);
    // A = A.T (transpot)

    let mut n = 0;
    let mut r = VecPoly::<K>::default();
    for i in 0..K
    {
        let mut tmp = [0u8; 32];
        let mut seed = [0u8; 32+2];
        seed[..32].copy_from_slice(&rc[..32]);
        seed[32..].copy_from_slice(&[n as u8, (n>>8) as u8]);
        prf(&seed, &mut tmp);
        r.poly[i] = CBD(&tmp, ETA1);
        n+=1;
    }
    let mut e1 = VecPoly::<K>::default();
    for i in 0..K
    {
        let mut tmp = [0u8; 32];
        let mut seed = [0u8; 32+2];
        seed[..32].copy_from_slice(&rc[..32]);
        seed[32..].copy_from_slice(&[n as u8, (n>>8) as u8]);
        prf(&seed, &mut tmp);
        e1.poly[i] = CBD(&tmp, ETA2);
        n+=1;
    }
    let mut tmp = [0u8; 32];
    let mut seed = [0u8; 32+2];
    seed[..32].copy_from_slice(&rc[..32]);
    seed[32..].copy_from_slice(&[n as u8, (n>>8) as u8]);
    prf(&seed, &mut tmp);
    let e2 = CBD(&tmp, ETA2);

    r = Ntt(&r);
    let mut u = inv_Ntt(&m_mult_v(&A, &Ntt(&r)));
    u = v_add(&u, &e1);

    let mut v = inv_ntt(&v_mult_v(&t, &r)); 
    add(&mut v, &e2);

    let mut tmp = Poly::default();
    decode(&mut tmp, m);
    tmp = decompress(&tmp, 1);
    
    add(&mut v, &tmp);


    let mut tmp = Compress(&u, Du as i16);
    Encode(&mut c, tmp, Du);
    let mut tmp = compress(&v, Dv as i16);
    encode(&mut c[CIPHERTEXTBYTES/2..], tmp);

    c
}


pub fn decryption(sk: &[u8], c: &[u8]) -> [u8; MESSAGEBYTES]
{
    let mut m = [0u8; MESSAGEBYTES];

    let mut tmp = VecPoly::<K>::default();
    Decode(&mut tmp, c, Du);
    let u = Decompress(&mut tmp, Du as i16);

    let mut tmp = Poly::default();
    decode(&mut tmp, c);//+Du *K*N/8);//, Du);
    let mut v = decompress(&mut tmp, Du as i16);

    let mut s = VecPoly::<K>::default();
    unpack_sk(sk, &mut s);


    let su = inv_ntt(&v_mult_v(&s, &Ntt(&u)));
    let su = neg(&su);
    add(&mut v, &su);

    let tmp = compress(&v, 1);
    encode(&mut m, tmp);

    // # Split ciphertext to vectors
    // index = self.du * self.k * self.R.n // 8
    // c2 = c[index:]
    
    // # Recover the vector u and convert to NTT form
    // u = self.M.decode(c, self.k, 1, l=self.du).decompress(self.du)
    // u.to_ntt()
    
    // # Recover the polynomial v
    // v = self.R.decode(c2, l=self.dv).decompress(self.dv)



    m
}