
use crate::params::*;
/* helper functions from the paper, bits and hints */

const Q_INV: i32 = 62209; // q^(-1) mod 2^16


/* for finite field element a with -2^{15}Q <= a <= Q*2^15,
 compute a' = a*2^{-16} (mod Q) such that -Q < a' < Q.
 https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
 here R = 2^16 */
pub fn montgomery(a: i32) -> i16
{
    let mut m = (a.wrapping_mul(Q_INV) as i16) as i32; 
    // a.wrapping_mul(b:T) returns (a * b) mod 2N, where N is the width of T in bits.
    m = (a - m * (Q as i32)) >> 16;
    assert!((m as i16) < Q as i16);
    m as i16
}

pub fn barrett(a: i16) -> i16 {
    let v = ((1u32 << 26) / Q as u32 + 1) as i32;
    let mut t = v * a as i32 + (1 << 25);
    t >>= 26;
    t *= Q as i32;
    a - t as i16
}


/* central reduction, returns r mod +- a */
pub fn cmod(r : i16, a: i16) -> i16
{
    let mut n = r % a;
    let mut t = a;
    if n > (t >> 1)
    {
        n -= t;
    }
    assert!(n <= (a+1)>>1);
    n 
}

pub fn mod_a(x: i16, a: i16) -> i16
{
    assert!((x%a + a) % a >= 0);
    return (x%a + a) % a;
}

pub fn _compress(x: i16, d: i16) -> i16
{
    let x = mod_a(x, Q as i16) as i32;
    let tmp: f32 = ((x<<(d as i32)) as f32)/(Q as f32);
    let tmp: i16 = tmp.round() as i16;
    return mod_a(tmp, 1<<d) as i16;
}

pub fn _decompress(x: i16, d: i16) -> i16
{
    let x = mod_a(x, Q as i16) as i32;
    let tmp: f32 = ((x as f32)*(Q as f32))/((1<<(d as i32)) as f32);
    return tmp.round() as i16 
}

