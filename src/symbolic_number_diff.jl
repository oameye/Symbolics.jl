# A complex-valued expression differentiated with respect to an ordinary symbolic
# independent variable remains one symbolic scalar. Non-holomorphic differentiation with
# respect to a complex independent variable is a separate calculus contract; this method
# only preserves the expression representation at the `Differential` application layer.
(D::Differential)(x::SymbolicNumber) = wrap(D(unwrap(x)))
