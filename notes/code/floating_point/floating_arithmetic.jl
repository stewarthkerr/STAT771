printstyled("\n********************************\n",color=:magenta)
printstyled("Finite Precision Arithmetic\n",color=:magenta)
printstyled("********************************\n\n",color=:magenta)

printstyled("Here, we add 10^3 + 2^(-6) in single precision.\n",color=:yellow)
printstyled("We expect the answer: "); printstyled("1000.000002\n",color=:blue)
printstyled("But, we actually get: "); printstyled("$(1f3 + 2f-6)\n",color=:red)
printstyled("It is as if, we did nothing.\n\n",color=:yellow)
