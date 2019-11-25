### Question 2    
    function listNormalizedFloatingPoint(base, digits, emin, emax)
        @assert emax >= emin
        @assert base >= 2

        #Initialize vector and put minimum normalized float inside
        float_vec = zeros(0)
        number = (float(base)^emin)
        append!(float_vec, number)

        #Loop through exponents
        for e = emin:emax
            #Loop through digits
            for d = digits:-1:1
                #Loop through base -- starting at 1 in so that we don't duplicated the previous number
                for b = 1:(base-1)
                    number = number + (float(base)^(e-d)*b)
                    append!(float_vec,number)
                end
            end
        end

        return float_vec

    end