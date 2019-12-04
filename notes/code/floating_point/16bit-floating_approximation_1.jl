printstyled("\n********************************\n",color=:magenta)
printstyled("IEEE 16 Bit Approximation near 1\n",color=:magenta)
printstyled("********************************\n\n",color=:magenta)

#Representation of 1
printstyled("First, we verify that the bitstring for 1 is ",color=:yellow)
    printstyled("0 01111 0000000000\n", color=:light_red)
    printstyled("The bitstring is: ")
    printstyled("$(x = bitstring(Float16(1)))\n\n", color=:green)

#Representation of next value
printstyled("We calculated next larger representable value after 1 to be ",color=:yellow)
    printstyled("1 + 2⁻¹⁰\n", color=:light_red)
    printstyled("Which has bitstring: ")
    printstyled("0 01111 0000000001\n", color=:light_red)
    printstyled("The bitstring for the next number is: ")
    printstyled("$(y = bitstring(nextfloat(Float16(1))))\n\n", color=:green)

#Floating Point Approximation (Lower half)
printstyled("What is the representation of any number from 1 to 1 + 2⁻¹⁰/2?\n", color=:yellow)
    printstyled("We claimed it should be ")
    printstyled("$x\n", color=:light_red)
    printstyled("For "); printstyled("$(z = 1 + 0.5*rand()*2^(-10))",color=:blue), printstyled(", the representation is ")
    printstyled("$(bitstring(Float16(z)  ))\n\n", color=:green)

#Floating point approximation (upper half)
printstyled("What is the representation of any number from 1 + 2⁻¹⁰/2 (not inclusive) to 1 + 2⁻¹⁰?\n", color=:yellow)
    printstyled("We claimed it should be ")
    printstyled("$y\n", color=:light_red)
    printstyled("For "); printstyled("$(z = 1 + 0.5*(1+rand())*2^(-10))",color=:blue), printstyled(", the representation is ")
    printstyled("$(bitstring(Float16(z)  ))\n\n", color=:green)
