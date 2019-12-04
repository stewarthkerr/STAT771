#IEEE 16 Bit Standard
printstyled("\n*********************************\n",color=:magenta)
printstyled("IEEE 16 Bit Approximation near 64\n",color=:magenta)
printstyled("*********************************\n\n",color=:magenta)

#Representation of 64
printstyled("First, we verify that the bitstring for 64 is ",color=:yellow)
    printstyled("0 10101 0000000000\n", color=:light_red)
    printstyled("The bitstring is: ")
    printstyled("$(x = bitstring(Float16(64)))\n\n", color=:green)

#Representation of next value
printstyled("We calculated next larger representable value after 64 to be ",color=:yellow)
    printstyled("64.0625\n", color=:light_red)
    printstyled("Which has bitstring: ")
    printstyled("0 10101 0000000001\n", color=:light_red)
    printstyled("The bitstring for the next number is: ")
    printstyled("$(y = bitstring(nextfloat(Float16(64))))\n\n", color=:green)

#Floating Point Approximation (Lower half)
printstyled("What is the representation of any number from 64 to 64 + 0.0625/2?\n", color=:yellow)
    printstyled("We claimed it should be ")
    printstyled("$x\n", color=:light_red)
    printstyled("For "); printstyled("$(z = 64 + 0.5*rand()*0.0625)",color=:blue), printstyled(", the representation is ")
    printstyled("$(bitstring(Float16(z)  ))\n\n", color=:green)

#Floating point approximation (upper half)
printstyled("What is the representation of any number from 64 + 0.0625/2 (not inclusive) to 64.0625?\n", color=:yellow)
    printstyled("We claimed it should be ")
    printstyled("$y\n", color=:light_red)
    printstyled("For "); printstyled("$(z = 64 + 0.5*(1+rand())*0.0625)",color=:blue), printstyled(", the representation is ")
    printstyled("$(bitstring(Float16(z)  ))\n\n", color=:green)

###################################################################
