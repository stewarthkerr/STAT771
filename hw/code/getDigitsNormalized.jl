### Question 2
    function getDigitsNormalized(decimalNum,base,digits)
        @assert base >= 2

        #Create a list that will hold the digits
        dig_list = zeros(Int, digits)

        #First, we find the best exponent
        exponent = trunc(Int,log(base,decimalNum))

        #Now, we find each of the digits for this exponent
        remainder = decimalNum
        for i = 0:(digits-1)
            #Find the value for the ith digit
            x = trunc(Int,remainder / (float(base)^(exponent-i)))
            dig_list[i+1] = x 

            #Calculate remainder for next iteration
            remainder = remainder % (float(base)^(exponent-i))
            print(i,' ',remainder,' ')
        end

        return (dig_list, exponent)
    end
