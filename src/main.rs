use bls12_381::*;
use ff::Field;
use group::{Curve, Group};
use rand::thread_rng;

fn spk_verify(
    w: G2Affine,
    a_prime: G1Affine,
    a_bar: G1Projective,
    d: G1Affine,
    sdk: (Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, G1Affine, G1Affine),
    revealed_msg: Scalar,
    generators: (G1Projective, G2Projective, G1Affine, G1Affine, G1Affine)
) -> bool {
    let (c, e_hat, r2_hat, r3_hat, s_hat, hidden_msg_hat, c1, c2) = sdk;
    let (p1, p2, h0, h_1, h_2) = generators;

    let c1_v =  (((a_bar - d) * c) + a_prime * e_hat + h0 * r2_hat).to_affine();
    let t = p1 + h_1 * revealed_msg;
    let c2_v = (t * c + d * (-r3_hat) + h0 * s_hat + h_2 * hidden_msg_hat).to_affine();

    if c1 != c1_v {return false};
    if c2 != c2_v {return false};
    if pairing(&(a_bar.to_affine()), &(p2.to_affine())) != pairing(&a_prime, &w)
    {
        return false
    }
    true
}


fn main() {
    let mut rng = thread_rng();
    // base points
    let p1 = G1Projective::generator();
    let p2 = G2Projective::generator();
    // issuer keys
    let x = Scalar::random(&mut rng);
    let w = (p2 * x).to_affine();
    // message generators
    let h0 = G1Projective::random(&mut rng).to_affine();
    let h_1 = G1Projective::random(&mut rng).to_affine();
    let h_2 = G1Projective::random(&mut rng).to_affine();
    // messages
    let m_1 = Scalar::random(&mut rng);
    let m_2 = Scalar::random(&mut rng);
    // signature
    let e = Scalar::random(&mut rng);
    let s = Scalar::random(&mut rng);
    let b = (p1 + h0 * s + h_1 * m_1 + h_2 * m_2).to_affine();
    let a = (b * (e + x).invert().unwrap()).to_affine();

    // Check regular signature
    assert_eq!(
        pairing(&b, &p2.to_affine()),
        pairing(&a, &(w + p2 * e).to_affine())
    );

    // Normal proof
    let r1 = Scalar::random(&mut rng);
    let r2 = Scalar::random(&mut rng);

    let r3 = r1.invert().unwrap();

    let e_tilde = Scalar::random(&mut rng);
    let s_tilde = Scalar::random(&mut rng);
    let r2_tilde = Scalar::random(&mut rng);
    let r3_tilde = Scalar::random(&mut rng);
    let m_2_tilde = Scalar::random(&mut rng);

    let a_prime = (a * r1).to_affine();
    let a_bar = a_prime * (-e) + b * r1;
    let d = (b * r1 + h0 * r2).to_affine();
    let s_prime = s + r2*r3;

    let c1_p = (a_prime * e_tilde + h0 * r2_tilde).to_affine();
    let c2_p = (d * (-r3_tilde) + h0 * s_tilde + h_2 * m_2_tilde).to_affine();

    let c = Scalar::random(&mut rng);

    let e_hat = e_tilde + c * e;
    let s_hat = s_tilde + c * s_prime;
    let r2_hat = r2_tilde + c * r2;
    let r3_hat = r3_tilde + c * r3;
    let m_2_hat = m_2_tilde + c * m_2;

    let sdk = (c, e_hat, r2_hat, r3_hat, s_hat, m_2_hat, c1_p, c2_p);
    let generators = (p1, p2, h0, h_1, h_2);
    let result = spk_verify(w, a_prime, a_bar, d, sdk, m_1, generators);

    assert_eq!(result, true);

    // Proof validated with w * y
    let y = Scalar::random(&mut rng);
    let w_f = (w * y).to_affine();

    let a_bar_f = a_prime * (-e) * y + b  * r1 * y;
    let d_f = (b * r1 * y + h0 * r2).to_affine();

    let c1_p_f = (a_prime * (e_tilde) + h0 * r2_tilde).to_affine();
    let c2_p_f = (d_f * (-r3_tilde) + h0 * s_tilde + h_2 * m_2_tilde).to_affine();

    let r3_hat_f = r3_tilde + c * r3 * y.invert().unwrap();
    let e_hat_f = e_tilde + c * y * e;
    let s_prime_f = s + r2 * r3 * y.invert().unwrap();
    let s_hat_f = s_tilde + c * s_prime_f;

    let sdk_f = (c, e_hat_f, r2_hat, r3_hat_f, s_hat_f, m_2_hat, c1_p_f, c2_p_f);
    let result_f = spk_verify(w_f, a_prime, a_bar_f, d_f, sdk_f, m_1, generators);

    assert_eq!(result_f, true);
}