printstyled("\n********************************\n",color=:magenta)
printstyled("ULP and Machine Epsilon in Julia\n",color=:magenta)
printstyled("********************************\n\n",color=:magenta)

#ULP
printstyled("For a number x, 1 ULP can be found using `eps(x)`\n\n",color=:yellow)
printstyled("Examples\tx\t\t\t\t\t\t1 ULP\n",color=:blue)
printstyled("--------\t-----\t\t\t\t\t------\n")
printstyled("Float16\t\t$(x=randn(Float16))\t\t\t\t\t$(eps(x))\n")
printstyled("Float32\t\t$(x=randn(Float32))\t\t\t\t$(eps(x))\n")
printstyled("Float64\t\t$(x=randn(Float64))\t\t$(eps(x))\n\n")

#Machine Epsilon
printstyled("Machine Epsilon corresponds to 1 ULP for 1.0\n\n",color=:yellow)
printstyled("Precision\tMachine Epsilon\n",color=:blue)
printstyled("Float16\t\t$(eps(Float16(1.0)))\n")
printstyled("Float32\t\t$(eps(1f0))\n")
printstyled("Float64\t\t$(eps())\n\n")
