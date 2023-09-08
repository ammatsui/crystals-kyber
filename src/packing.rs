use crate::{params::*, poly::*,};
use bitvec::prelude::*;


/* the same as pack/unpack from dilithium */
/* 256 coeff -> each into l bits -> glue together and get into 32*l bytes */
pub fn encode(r: &mut [u8], a: &Poly, l: usize)
{
    for (slot, c) in r[..32*l]
  .view_bits_mut::<Lsb0>()
  .chunks_mut(l)
  .zip(a.coeff.iter().copied())
  {
    slot.store_le(c as u16); // because the coeffs are supposed to be positive
    assert_eq!(slot.load_le::<u16>(), c as u16); 
  }
} 


/* array of 32*l bytes -> 256*l bits -> l bits into one coeff -> 256 coeff of a polynomial */
pub fn decode(r: &mut Poly, a: &[u8], l: usize)
{
    let bits = a[..32*l].view_bits::<Lsb0>();
    for i in 0..N 
    {
        r.coeff[i] = bits[l*i..l*(i+1)].load_le::<u16>() as i16;
    }
}


pub fn Encode(r: &mut [u8], a: &VecPoly<K>, d: usize) 
{
    for i in 0..K 
    {
        encode(&mut r[ i * d*N/8..], &a.poly[i], d);
    }
}


/* decode */
pub fn Decode(r: &mut VecPoly<K>, a: &[u8], d: usize) {
    for i in 0..K 
    {
        decode( &mut r.poly[i], &a[i * d*N/8..], d);
    }
}


/* main packers */
/* secret key sk = (s)*/
pub fn pack_sk(sk: &mut [u8], s: &VecPoly<K>)
{
    for i in 0..K 
    {
        encode(&mut sk[ i * PACKED_KEYS*N/8..], &s.poly[i], PACKED_KEYS);
    }
}


pub fn unpack_sk( sk: &[u8], s: &mut VecPoly<K>)
{
    for i in 0..K 
    {
        decode( &mut s.poly[i], &sk[i * PACKED_KEYS*N/8..], PACKED_KEYS);
        s.poly[i].ntt = true;
    }
  }
  


/* public key pk = (t, rho) */
pub fn pack_pk(pk: &mut [u8], t: &VecPoly<K>, rho: &[u8])
{
    pk[..SEEDBYTES].copy_from_slice(&rho[..SEEDBYTES]);
    for i in 0..K 
    {
        encode(&mut pk[SEEDBYTES + i * PACKED_KEYS*N/8..], &t.poly[i], PACKED_KEYS);
    }
}

/* unpack public key pk = (t, rho) */
pub fn unpack_pk(pk: &[u8], t: &mut VecPoly<K>, rho: &mut [u8])
{
    rho[..SEEDBYTES].copy_from_slice(&pk[..SEEDBYTES]);
    for i in 0..K 
    {
        decode( &mut t.poly[i], &pk[SEEDBYTES + i * PACKED_KEYS*N/8..], PACKED_KEYS);
        t.poly[i].ntt = true;
    }
}