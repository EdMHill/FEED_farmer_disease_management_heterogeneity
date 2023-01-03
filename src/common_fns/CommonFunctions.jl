#=
Purpose:
File to house other commonly used functions

Alogorithm list:
 - Calculation for 1-exp(x)
=#

#-------------------------------------------------------------------------------
### Code for calculation of "1 - exp(x)"
#-------------------------------------------------------------------------------

"""
    oneMinusExp(x::Float64)

Return solution to "1 - exp(x)" using algorithm described at http://www.johndcook.com/blog/cpp_expm1/. \n
Done as if x is very small, directly computing 1 - exp(x) can be inaccurate.

Location: CommonFunctions.jl
"""
function oneMinusExp(x::Float64)

    if x == 0.0
        return 0.0
    elseif abs(x) < 1e-5  #If x is of small enough order, use taylor expansion of exp(x).
                        # Applying to 1 - exp(x), keep first couple of remaining terms -(x + (x^2/2))
        return -(x + 0.5*x*x)
    else
        return -(exp(x) - 1.0)
    end

end
