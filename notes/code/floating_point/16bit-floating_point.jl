
#IEEE 16 Bit Standard
printstyled("\n********************\n",color=:magenta)
printstyled("IEEE 16 Bit Standard\n",color=:magenta)
printstyled("********************\n\n",color=:magenta)

#Number of Digits
printstyled("Here, we verify that 16-bit numbers are stored with 16 bits.\n",color=:yellow)
printstyled("First, we generate a random 16-bit number: ");
    printstyled("$(x=randn(Float16))\n", color=:green)
printstyled("Second, we ask Julia for its binary representation: ")
    printstyled("$(bl = bitstring(x))\n", color=:green)
printstyled("Finally, we check the length of the string: ")
    printstyled("$(length(bl))\n\n", color=:green)

#Representing Negative vs Positive Numbers
printstyled("The first bit of -1.0 and 1.0 are different.\n",color=:yellow)
printstyled("-1 is: ")
    printstyled("$(bitstring(Float16(-1.0)))\n",color=:green)
printstyled("1 is: ")
    printstyled("$(bitstring(Float16(1.0)))\n", color=:green)
printstyled("The leading digit is 1 for -1, and 0 for 1.\n\n",color=:yellow)

#Largest Positive 16-bit floating point number
printstyled("We computed the largest number to be 65504.\nWe now verify this with Julia.\n",color=:yellow)
    printstyled("Largest 16-bit float in Julia: ")
    printstyled("$(x = prevfloat(Inf16))\n",color=:green)
    printstyled("Is this equal to 65504? ")
    printstyled("$(65504 == x)",color=:green)

#Smallest positive 16-bit normalized floating point number
printstyled("We compute the smallest positive normalized 16-bit number.\n",color=:yellow)
    printstyled("The binary representation is: ")
    printstyled("$(x = "0000010000000000")\n",color=:green)
    printstyled("This is equal to: ")
    printstyled("$(y = bitstring(Float16(2^(-14))))\n\n",color=:green)

#Smallest 16-bit subnormal floating point number
printstyled("We compute the smallest positive subnormal 16-bit number.\n",color=:yellow)
    printstyled("We believe the number is: ")
    printstyled("$(x = Float16(2^(-24)))\n", color=:green)
    printstyled("The actual number is: ")
    printstyled("$(y = nextfloat(Float16(0)))\n",color=:green)
    printstyled("Are these equal? ")
    printstyled("$(x==y)\n\n", color=:green)
