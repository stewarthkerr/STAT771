#IEEE 32 Bit Standard
printstyled("\n********************\n",color=:magenta)
printstyled("IEEE 32 Bit Standard\n",color=:magenta)
printstyled("********************\n\n",color=:magenta)

#Largest Positive 32-bit floating point number
printstyled("We computed the largest number to be 2¹²⁸( 1 - 2⁻²⁴) .\nWe now verify this with Julia.\n",color=:yellow)
    printstyled("Largest 32-bit float in Julia: ")
    printstyled("$(x = prevfloat(Inf32))\n",color=:green)
    printstyled("Is this equal to 2¹²⁸( 1 - 2⁻²⁴)? ")
    printstyled("$(2.0^128*(1 - 2.0^(-24)) == x)\n\n",color=:green)

#Smallest positive 16-bit normalized floating point number
printstyled("We compute the smallest positive normalized 32-bit number.\n",color=:yellow)
    printstyled("The binary representation is: ")
    printstyled("$(x = "00000000100000000000000000000000")\n",color=:green)
    printstyled("This is equal to: ")
    printstyled("$(y = bitstring(Float32(2^(-126))))\n\n",color=:green)

#Smallest 16-bit subnormal floating point number
printstyled("We compute the smallest positive subnormal 32-bit number.\n",color=:yellow)
    printstyled("We believe the number is: ")
    printstyled("$(x = Float32(2^(-149)))\n", color=:green)
    printstyled("The actual number is: ")
    printstyled("$(y = nextfloat(Float32(0)))\n",color=:green)
    printstyled("Are these equal? ")
    printstyled("$(x==y)\n\n", color=:green)
