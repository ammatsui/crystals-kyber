use crate::{params::*, utils::*};
/* ntt library: forward, inverse,  multiplication */

/* roots of unity in order for the forward ntt */
pub const ZETAS: [i16; 128] = [2285, 2571, 2970, 1812, 1493, 1422, 287, 202, 3158, 622, 1577, 182, 962, 2127, 1855, 1468, 
573, 2004, 264, 383, 2500, 1458, 1727, 3199, 2648, 1017, 732, 608, 1787, 411, 3124, 1758, 
1223, 652, 2777, 1015, 2036, 1491, 3047, 1785, 516, 3321, 3009, 2663, 1711, 2167, 126, 1469, 
2476, 3239, 3058, 830, 107, 1908, 3082, 2378, 2931, 961, 1821, 2604, 448, 2264, 677, 2054, 
2226, 430, 555, 843, 2078, 871, 1550, 105, 422, 587, 177, 3094, 3038, 2869, 1574, 1653, 3083, 
778, 1159, 3182, 2552, 1483, 2727, 1119, 1739, 644, 2457, 349, 418, 329, 3173, 3254, 817, 
1097, 603, 610, 1322, 2044, 1864, 384, 2114, 3193, 1218, 1994, 2455, 220, 2142, 1670, 2144, 
1799, 2051, 794, 1819, 2475, 2459, 478, 3221, 3021, 996, 991, 958, 1869, 1522, 1628];


pub fn ntt_(a: &mut [i16])
/* arg: array of coefficients 
output in bit-reversed order */
{
    let mut j;
    let mut k = 1;
    let mut len = N/2; //128
    let (mut t, mut zeta);
  
    while len >= 2 {
      let mut i = 0;
      while i < N {
        zeta = ZETAS[k] as i32;
        k += 1;
        j = i;
        while j < (i + len) {
          t = montgomery(zeta * (a[j + len] as i32));
          a[j + len] = a[j] - t;
          a[j] += t;
          j += 1;
        }
        i = j + len;
      }
      len >>= 1;
    }

    for i in 0..N 
    {
      assert!(a[i].abs() < Q as i16);
    }
  }


/* inverse ntt (and mult by Montgomery factor r = 2^32) */
pub fn inv_ntt_(a: &mut [i16])
{
    let mut j;
    let mut k = N/2 - 1;
    let mut len = 2;
    let (mut t, mut zeta);
    const F: i32 = 1441; // mont^2/256
  
    while len <= N/2 {
      let mut i = 0;
      while i < N {
        zeta = -ZETAS[k] as i32;
        k -= 1;
        j = i;
        while j < (i + len) {
          t = a[j];
          a[j] = barrett(t + a[j + len]);
          a[j + len] = a[j+len] - t;//t - a[j + len];
          a[j + len] = montgomery(zeta * a[j + len] as i32);
          j += 1
        }
        i = j + len;
      }
      len <<= 1;
    }
    for j in 0..N {
      a[j] = montgomery(F * a[j] as i32);
    }
  }
  