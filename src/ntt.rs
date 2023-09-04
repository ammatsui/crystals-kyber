use crate::{params::*, utils::*};


pub const ZETAS: [i16; 128] = [
    -1044, -758, -359, -1517, 1493, 1422, 287, 202, -171, 622, 1577, 182, 962, -1202, -1474, 1468,
    573, -1325, 264, 383, -829, 1458, -1602, -130, -681, 1017, 732, 608, -1542, 411, -205, -1571,
    1223, 652, -552, 1015, -1293, 1491, -282, -1544, 516, -8, -320, -666, -1618, -1162, 126, 1469,
    -853, -90, -271, 830, 107, -1421, -247, -951, -398, 961, -1508, -725, 448, -1065, 677, -1275,
    -1103, 430, 555, 843, -1251, 871, 1550, 105, 422, 587, 177, -235, -291, -460, 1574, 1653, -246,
    778, 1159, -147, -777, 1483, -602, 1119, -1590, 644, -872, 349, 418, 329, -156, -75, 817, 1097,
    603, 610, 1322, -1285, -1465, 384, -1215, -136, 1218, -1335, -874, 220, -1187, -1659, -1185,
    -1530, -1278, 794, -1510, -854, -870, 478, -108, -308, 996, 991, 958, -1460, 1522, 1628,
];

/// Name:  fqmul
///
/// Description: Multiplication followed by Montgomery reduction
///
/// Arguments:   - i16 a: first factor
///  - i16 b: second factor
///
/// Returns 16-bit integer congruent to a*b*R^{-1} mod q
pub fn fqmul(a: i16, b: i16) -> i16 {
    montgomery(a as i32 * b as i32)
}

/// Name:  ntt
///
/// Description: Inplace number-theoretic transform (NTT) in Rq
///  input is in standard order, output is in bitreversed order
///
/// Arguments:   - i16 r[256]: input/output vector of elements of Zq
pub fn ntt_(r: &mut [i16]) {
    let mut j;
    let mut k = 1usize;
    let mut len = 128;
    let (mut t, mut zeta);

    while len >= 2 {
        let mut start = 0;
        while start < 256 {
            zeta = ZETAS[k];
            k += 1;
            j = start;
            while j < (start + len) {
                t = fqmul(zeta, r[j + len]);
                r[j + len] = r[j] - t;
                r[j] += t;
                j += 1;
            }
            start = j + len;
        }
        len >>= 1;
    }

}

/// Name:  invntt
///
/// Description: Inplace inverse number-theoretic transform in Rq
///  input is in bitreversed order, output is in standard order
///
/// Arguments:   - i16 r[256]: input/output vector of elements of Zq
pub fn inv_ntt_(r: &mut [i16]) {
    let mut j;
    let mut k = 127usize;
    let mut len = 2;
    let (mut t, mut zeta);
    const F: i16 = 1441; // mont^2/128
    while len <= 128 {
        let mut start = 0;
        while start < 256 {
            zeta = ZETAS[k];
            k -= 1;
            j = start;
            while j < (start + len) {
                t = r[j];
                r[j] = barrett(t + r[j + len]);
                r[j + len] = r[j + len] - t;
                r[j + len] = fqmul(zeta, r[j + len]);
                j += 1
            }
            start = j + len;
        }
        len <<= 1;
    }
    for j in 0..256 {
        r[j] = fqmul(r[j], F);
    }
}












// use crate::{params::*, utils::*};
// /* ntt library: forward, inverse,  multiplication */

// /* roots of unity in order for the forward ntt */
// // pub const ZETAS: [i16; 128] = [2285, 2571, 2970, 1812, 1493, 1422, 287, 202, 3158, 622, 1577, 182, 962, 2127, 1855, 1468, 
// // 573, 2004, 264, 383, 2500, 1458, 1727, 3199, 2648, 1017, 732, 608, 1787, 411, 3124, 1758, 
// // 1223, 652, 2777, 1015, 2036, 1491, 3047, 1785, 516, 3321, 3009, 2663, 1711, 2167, 126, 1469, 
// // 2476, 3239, 3058, 830, 107, 1908, 3082, 2378, 2931, 961, 1821, 2604, 448, 2264, 677, 2054, 
// // 2226, 430, 555, 843, 2078, 871, 1550, 105, 422, 587, 177, 3094, 3038, 2869, 1574, 1653, 3083, 
// // 778, 1159, 3182, 2552, 1483, 2727, 1119, 1739, 644, 2457, 349, 418, 329, 3173, 3254, 817, 
// // 1097, 603, 610, 1322, 2044, 1864, 384, 2114, 3193, 1218, 1994, 2455, 220, 2142, 1670, 2144, 
// // 1799, 2051, 794, 1819, 2475, 2459, 478, 3221, 3021, 996, 991, 958, 1869, 1522, 1628];
// pub const ZETAS: [i16; 128] = [
//   -1044, -758, -359, -1517, 1493, 1422, 287, 202, -171, 622, 1577, 182, 962, -1202, -1474, 1468,
//   573, -1325, 264, 383, -829, 1458, -1602, -130, -681, 1017, 732, 608, -1542, 411, -205, -1571,
//   1223, 652, -552, 1015, -1293, 1491, -282, -1544, 516, -8, -320, -666, -1618, -1162, 126, 1469,
//   -853, -90, -271, 830, 107, -1421, -247, -951, -398, 961, -1508, -725, 448, -1065, 677, -1275,
//   -1103, 430, 555, 843, -1251, 871, 1550, 105, 422, 587, 177, -235, -291, -460, 1574, 1653, -246,
//   778, 1159, -147, -777, 1483, -602, 1119, -1590, 644, -872, 349, 418, 329, -156, -75, 817, 1097,
//   603, 610, 1322, -1285, -1465, 384, -1215, -136, 1218, -1335, -874, 220, -1187, -1659, -1185,
//   -1530, -1278, 794, -1510, -854, -870, 478, -108, -308, 996, 991, 958, -1460, 1522, 1628,
// ];


// pub fn ntt_(a: &mut [i16])
// /* arg: array of coefficients 
// output in bit-reversed order */
// {
//     println!("{:?}", a);
//     let mut j;
//     let mut k = 1;
//     let mut len = N/2; //128
//     let (mut t, mut zeta);
  
//     while len >= 2 {
//       let mut i = 0;
//       while i < N {
//         zeta = ZETAS[k] as i32;
//         k += 1;
//         j = i;
//         while j < (i + len) {
//           t = montgomery(zeta * (a[j + len] as i32));
//           a[j + len] = a[j] - t;
//           a[j] += t;
//           j += 1;
//         }
//         i = j + len;
//       }
//       len >>= 1;
//     }

//     println!("*************");
//     println!("{:?}", a);
//     for i in 0..N 
//     {
//       assert!(a[i].abs() <= Q as i16);
//     }
//   }


// /* inverse ntt (and mult by Montgomery factor r = 2^32) */
// pub fn inv_ntt_(a: &mut [i16])
// {
//     let mut j;
//     let mut k = N/2 - 1;
//     let mut len = 2;
//     let (mut t, mut zeta);
//     const F: i32 = 1441; // mont^2/256
  
//     while len <= N/2 {
//       let mut i = 0;
//       while i < N {
//         zeta = -ZETAS[k] as i32;
//         k -= 1;
//         j = i;
//         while j < (i + len) {
//           t = a[j];
//           a[j] = barrett(t + a[j + len]);
//           a[j + len] = a[j+len] - t;//t - a[j + len];
//           a[j + len] = montgomery(zeta * a[j + len] as i32);
//           j += 1
//         }
//         i = j + len;
//       }
//       len <<= 1;
//     }
//     for j in 0..N {
//       a[j] = montgomery(F * a[j] as i32);
//     }
//   }
  