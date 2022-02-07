pragma circom 2.0.2;

include "../../circuits/elliptic.circom";

component main {public [a, b]} = EllipticCurveAddUnequal4Reg(96,
							     54880396502181392957329877675,
							     31935979117156477062286671870,
							     20826981314825584179608359615,
							     8047903782086192180586325942);
