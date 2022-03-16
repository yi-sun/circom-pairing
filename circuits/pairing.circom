pragma circom 2.0.2;

include "bigint.circom";
include "fp.circom";
include "fp12.circom";
include "../node_modules/circomlib/circuits/comparators.circom";


// NOT USABLE YET. pretends point coordinates are just field elements for now, 
// because we don't have nonnative field arithmetic yet
// copied formulas from python/pairing: 
// def linefunc(P1, P2, T):
//     assert P1 and P2 and T # No points-at-infinity allowed, sorry
//     x1, y1 = P1
//     x2, y2 = P2
//     xt, yt = T
//     if x1 != x2:
//         m = (y2 - y1) / (x2 - x1)
//         return m * (xt - x1) - (yt - y1)
//     elif y1 == y2:
//         m = 3 * x1**2 / (2 * y1)
//         return m * (xt - x1) - (yt - y1)
//     else:
//         return xt - x1
template Linefunc(n, k) {
    signal input p1[2][k];
    signal input p2[2][k];
    signal input t[2][k];
    signal input p[k];
    signal output out[k];

    component EqualX = BigIsEqual(k);
    component EqualY = BigIsEqual(k);
    for (var i = 0; i < k; i ++) {
        EqualX.a[i] <== p1[0][i];
        EqualX.b[i] <== p2[0][i];
        EqualY.a[i] <== p1[1][i];
        EqualY.b[i] <== p2[1][i];
    }
    signal unequalX; // x1 != x2
    signal equalY; // y1 == y2
    signal equalYequalX; // y1 == y2 and x1 == x2
    signal unequalYequalX; // y1 != y2 and x1 == x2

    unequalX <== 1 - EqualX.out;
    equalY <== EqualY.out;
    equalYequalX <== equalY - unequalX * equalY;
    unequalYequalX <== 1 - unequalX - equalYequalX;

    component Y2Y1 = BigSubModP(n, k); // to compute y2-y1
    component X2X1 = BigSubModP(n, k); // to compute x2-x1
    component YTY1 = BigSubModP(n, k); // to compute yt-y1
    component XTX1 = BigSubModP(n, k); // to compute xt-x1

    for (var i = 0; i < k; i ++) {
        Y2Y1.a[i] <== p2[1][i];
        Y2Y1.b[i] <== p1[1][i];
        X2X1.a[i] <== p2[0][i];
        X2X1.b[i] <== p1[0][i];
        YTY1.a[i] <== t[1][i];
        YTY1.b[i] <== p1[1][i];
        XTX1.a[i] <== t[0][i];
        XTX1.b[i] <== p1[0][i];

        Y2Y1.p[i] <== p[i];
        X2X1.p[i] <== p[i];
        YTY1.p[i] <== p[i];
        XTX1.p[i] <== p[i];
    }
    signal y2y1[k];
    signal x2x1[k];
    signal yty1[k];
    signal xtx1[k];

    for (var i = 0; i < k; i ++) {
        y2y1[i] <== Y2Y1.out[i];
        x2x1[i] <== X2X1.out[i];
        yty1[i] <== YTY1.out[i];
        xtx1[i] <== XTX1.out[i];
    }
    component modInv1 = BigModInv(n, k); // to invert x2-x1
    component modInv2 = BigModInv(n, k); // to invert 2y1
    component DoubleY1 = BigAddModP(n, k); // to compute 2y1

    for (var i = 0; i < k; i ++) {
        DoubleY1.a[i] <== p1[1][i];
        DoubleY1.b[i] <== p1[1][i];
        DoubleY1.p[i] <== p[i];
    }
    signal doubleY1[k];
    for (var i = 0; i < k; i ++) {
        doubleY1[i] <== DoubleY1.out[i];
    }

    // if x2 = x1 we need nonzero residue, else we get BigModInv assertion error
    modInv1.in[0] <== x2x1[0] + (1 - unequalX);
    modInv1.p[0] <== p[0];
    modInv2.in[0] <== doubleY1[0];
    modInv2.p[0] <== p[0];

    for (var i = 1; i < k; i ++) {
        modInv1.in[i] <== x2x1[i];
        modInv1.p[i] <== p[i];
        modInv2.in[i] <== doubleY1[i];
        modInv2.p[i] <== p[i];
    }
    signal x2x1inv[k];
    signal doubleY1inv[k];
    for (var i = 0; i < k; i ++) {
        x2x1inv[i] <== modInv1.out[i];
        doubleY1inv[i] <== modInv2.out[i];
    }
    component SquareX1 = BigMultModP(n, k); // to compute x1*x1
    for (var i = 0; i < k; i ++) {
        SquareX1.a[i] <== p1[0][i];
        SquareX1.b[i] <== p1[0][i];
        SquareX1.p[i] <== p[i];
    }
    signal squareX1[k];
    for (var i = 0; i < k; i ++) {
        squareX1[i] <== SquareX1.out[i];
    }
    component TwoSquareX1 = BigAddModP(n, k); // to compute 2x1*x1
    for (var i = 0; i < k; i ++) {
        TwoSquareX1.a[i] <== squareX1[i];
        TwoSquareX1.b[i] <== squareX1[i];
        TwoSquareX1.p[i] <== p[i];
    }
    signal twoSquareX1[k];
    for (var i = 0; i < k; i ++) {
        twoSquareX1[i] <== TwoSquareX1.out[i];
    }
    component ThreeSquareX1 = BigAddModP(n, k); // to compute 3x1*x1
    for (var i = 0; i < k; i ++) {
        ThreeSquareX1.a[i] <== twoSquareX1[i];
        ThreeSquareX1.b[i] <== squareX1[i];
        ThreeSquareX1.p[i] <== p[i];
    }
    signal threeSquareX1[k];
    for (var i = 0; i < k; i ++) {
        threeSquareX1[i] <== ThreeSquareX1.out[i];
    }
    component Mult1 = BigMultModP(n, k); // to compute (y2 - y1) / (x2 - x1)
    component Mult2 = BigMultModP(n, k); // to compute 3 * x1**2 / (2 * y1)
    for (var i = 0; i < k; i ++) {
        Mult1.a[i] <== y2y1[i];
        Mult1.b[i] <== x2x1inv[i];
        Mult1.p[i] <== p[i];

        Mult2.a[i] <== threeSquareX1[i];
        Mult2.b[i] <== doubleY1inv[i];
        Mult2.p[i] <== p[i];
    }

    signal mult1[k];
    signal mult2[k];
    for (var i = 0; i < k; i ++) {
        mult1[i] <== Mult1.out[i];
        mult2[i] <== Mult2.out[i];
    }
    signal m1[k];
    signal m[k];
    for (var i = 0; i < k; i ++) {
        m1[i] <== unequalX * mult1[i];
        m[i] <== m1[i] + equalY * mult2[i];
    }

    component MX = BigMultModP(n, k); // to compute m * (xt - x1)
    for (var i = 0; i < k; i ++) {
        MX.a[i] <== m[i];
        MX.b[i] <== xtx1[i];
        MX.p[i] <== p[i];
    }

    signal mx[k];
    for (var i = 0; i < k; i ++) {
        mx[i] <== MX.out[i];
    }
    component MXY = BigSubModP(n, k); // to compute m * (xt - x1) - (yt - y1)
    for (var i = 0; i < k; i ++) {
        MXY.a[i] <== mx[i];
        MXY.b[i] <== yty1[i];
        MXY.p[i] <== p[i];
    }
    signal mxy[k];
    for (var i = 0; i < k; i ++) {
        mxy[i] <== MXY.out[i];
    }

    signal result[k];
    for (var i = 0; i < k; i ++) { // use m * (xt - x1) - (yt - y1) or xt - x1
        result[i] <== mxy[i] - unequalYequalX * mxy[i];
        out[i] <== result[i] + unequalYequalX * xtx1[i];
    }
}
