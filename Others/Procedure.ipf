#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function p_net()

	wave BetheCoefs
	
	If (BetheCoefs[4] > 0.5)
	BetheCoefs[3] = BetheCoefs[2] * 2 * (1 - BetheCoefs[4])
	Endif
	
	If (BetheCoefs[4] <= 0.5)
	BetheCoefs[3] = BetheCoefs[2] * 2 * BetheCoefs[4]
	Endif
	
	return BetheCoefs[2]
	
end

Function p_A()

	wave BetheCoefs

	If (BetheCoefs[4] == 0)	
	BetheCoefs[5] = 0
	Endif
	
	If (BetheCoefs[4] != 0)	
	BetheCoefs[5] = BetheCoefs[3] / (2 * BetheCoefs[4])
	Endif
		
	return BetheCoefs[5]

end

Function p_B()

	wave BetheCoefs

	If (BetheCoefs[4] == 1)	
	BetheCoefs[6] = 0
	Endif
	
	If (BetheCoefs[4] != 1)	
	BetheCoefs[6] = BetheCoefs[3] / (2 * (1 - BetheCoefs[4]))
	Endif
	
	return BetheCoefs[5]

end

Function R_fraction()

	wave BetheCoefs
	
	BetheCoefs[8] = (2 * (1 - BetheCoefs[7]) * BetheCoefs[4] + 2 * BetheCoefs[7] * BetheCoefs[4] * 2) / (2 * (1 - BetheCoefs[7]) * BetheCoefs[4] + 2 * BetheCoefs[7] * BetheCoefs[4] * 2 + 8 * (1 - BetheCoefs[4]))

End

Function NumberDensity()

	wave BetheCoefs
	
//	BetheCoefs[20] = (1 - BetheCoefs[7]) * BetheCoefs[4] * BetheCoefs[19] / ((1 - BetheCoefs[7]) * BetheCoefs[4] * BetheCoefs[17] + 2 * BetheCoefs[7] * BetheCoefs[4] * BetheCoefs[17] / 2 + (1 - BetheCoefs[7]) * BetheCoefs[18])
//	BetheCoefs[21] = (1 - BetheCoefs[4]) * BetheCoefs[19] / ((1 - BetheCoefs[7]) * BetheCoefs[4] * BetheCoefs[17] + 2 * BetheCoefs[7] * BetheCoefs[4] * BetheCoefs[17] / 2 + (1 - BetheCoefs[7]) * BetheCoefs[18])
//	BetheCoefs[22] = 2 * BetheCoefs[7] * BetheCoefs[4] * BetheCoefs[19] / ((1 - BetheCoefs[7]) * BetheCoefs[4] * BetheCoefs[17] + 2 * BetheCoefs[7] * BetheCoefs[4] * BetheCoefs[17] / 2 + (1 - BetheCoefs[7]) * BetheCoefs[18])
	
	BetheCoefs[20] = BetheCoefs[19] * BetheCoefs[24] * BetheCoefs[4] * (1 - BetheCoefs[7]) / (BetheCoefs[17] * BetheCoefs[24] * BetheCoefs[4] + BetheCoefs[18] * BetheCoefs[23] * (1 - BetheCoefs[4]))
	BetheCoefs[21] = BetheCoefs[19] * BetheCoefs[23] * (1 - BetheCoefs[4]) / (BetheCoefs[17] * BetheCoefs[24] * BetheCoefs[4] + BetheCoefs[18] * BetheCoefs[23] * (1 - BetheCoefs[4]))
	BetheCoefs[22] = BetheCoefs[19] * BetheCoefs[24] * BetheCoefs[4] * 2 * BetheCoefs[7] / (BetheCoefs[17] * BetheCoefs[24] * BetheCoefs[4] + BetheCoefs[18] * BetheCoefs[23] * (1 - BetheCoefs[4]))
	BetheCoefs[25] = BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22]

//	BetheCoefs[25] = BetheCoefs[19] / ((1 - BetheCoefs[7]) * BetheCoefs[4] * BetheCoefs[17] + 2 * BetheCoefs[7] * BetheCoefs[4] * BetheCoefs[17] / 2 + (1 - BetheCoefs[4]) * BetheCoefs[18])

End

Function Bethe()

	variable X_A8, X_A7, X_A6, X_A5, X_A4, X_A3, X_A2, X_A1, X_C4, X_C3, X_C2, X_C1, X_B8, X_B7, X_B6, X_B5, X_B4, X_B3, X_B2, X_B1
	variable X_A0, X_B0, X_C0
	variable N_A8, N_A7, N_A6, N_A5, N_A4, N_A3, N_A2, N_A1, N_C4, N_C3, N_C2, N_C1, N_B8, N_B7, N_B6, N_B5, N_B4, N_B3, N_B2, N_B1
	variable N_A0, N_B0, N_C0
	variable R_A, R_B, R_C
	wave BetheCoefs
	
	p_net()
	p_A()
	p_B()
	R_fraction()
	NumberDensity()
	
	BetheCoefs[9] = BetheCalc()
	
	If (BetheCoefs[5] == 1 && BetheCoefs[6] == 1)
	BetheCoefs[9] = 0
	EndIf
	
	If (BetheCoefs[9] < 0 && BetheCoefs[9] != -1)
	BetheCoefs[9] = 1
	Endif
	
	BetheCoefs[10] = 1 - BetheCoefs[6] + BetheCoefs[6] * ((1 - BetheCoefs[7]) * BetheCoefs[9]^(BetheCoefs[23] - 1) + BetheCoefs[7] * BetheCoefs[9]^(BetheCoefs[23] / 2 - 1))
	
	If (BetheCoefs[23] == 8)
	X_A8 = (1 - BetheCoefs[9])^8
	X_A7 = 8 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^7
	X_A6 = 28 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^6
	X_A5 = 56 * BetheCoefs[9]^3 * (1 - BetheCoefs[9])^5
	X_A4 = 70 * BetheCoefs[9]^4 * (1 - BetheCoefs[9])^4
	X_A3 = 56 * BetheCoefs[9]^5 * (1 - BetheCoefs[9])^3	
	X_A2 = 28 * BetheCoefs[9]^6 * (1 - BetheCoefs[9])^2
	X_A1 = 8 * BetheCoefs[9]^7 * (1 - BetheCoefs[9])^1
	Endif
	
	If (BetheCoefs[23] == 7)
	X_A8 = 0
	X_A7 = 1 * (1 - BetheCoefs[9])^7
	X_A6 = 7 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^6
	X_A5 = 21 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^5
	X_A4 = 35 * BetheCoefs[9]^3 * (1 - BetheCoefs[9])^4
	X_A3 = 35 * BetheCoefs[9]^4 * (1 - BetheCoefs[9])^3	
	X_A2 = 21 * BetheCoefs[9]^5 * (1 - BetheCoefs[9])^2
	X_A1 = 7 * BetheCoefs[9]^6 * (1 - BetheCoefs[9])^1
	Endif

	If (BetheCoefs[23] == 6)
	X_A8 = 0
	X_A7 = 0
	X_A6 = 1 * (1 - BetheCoefs[9])^6
	X_A5 = 6 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^5
	X_A4 = 15 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^4
	X_A3 = 20 * BetheCoefs[9]^3 * (1 - BetheCoefs[9])^3	
	X_A2 = 15 * BetheCoefs[9]^4 * (1 - BetheCoefs[9])^2
	X_A1 = 6 * BetheCoefs[9]^5 * (1 - BetheCoefs[9])^1
	Endif

	If (BetheCoefs[23] == 5)
	X_A8 = 0
	X_A7 = 0
	X_A6 = 0
	X_A5 = 1 * (1 - BetheCoefs[9])^5
	X_A4 = 5 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^4
	X_A3 = 10 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^3	
	X_A2 = 10 * BetheCoefs[9]^3 * (1 - BetheCoefs[9])^2	
	X_A1 = 5 * BetheCoefs[9]^4 * (1 - BetheCoefs[9])^1
	Endif

	If (BetheCoefs[23] == 4)
	X_A8 = 0
	X_A7 = 0
	X_A6 = 0
	X_A5 = 0
	X_A4 = 1 * (1 - BetheCoefs[9])^4
	X_A3 = 4 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^3	
	X_A2 = 6 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^2	
	X_A1 = 4 * BetheCoefs[9]^3 * (1 - BetheCoefs[9])^1	
	Endif
	
	If (BetheCoefs[23] == 3)
	X_A8 = 0
	X_A7 = 0
	X_A6 = 0
	X_A5 = 0
	X_A4 = 0
	X_A3 = 1 * (1 - BetheCoefs[9])^3	
	X_A2 = 3 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^2	
	X_A1 = 3 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^1
	Endif
	
	If (BetheCoefs[23] == 2)
	X_A8 = 0
	X_A7 = 0
	X_A6 = 0
	X_A5 = 0
	X_A4 = 0
	X_A3 = 0
	X_A2 = 1 * (1 - BetheCoefs[9])^2
	X_A1 = 2 * BetheCoefs[9]^1 * (1 - BetheCoefs[9])^1
	Endif
	
    If (BetheCoefs[23] == 8)
	X_C4 = (1 - BetheCoefs[9])^4
	X_C3 = 4 * BetheCoefs[9] * (1 - BetheCoefs[9])^3
	X_C2 = 6 * BetheCoefs[9]^2 * (1 - BetheCoefs[9])^2
	X_C1 = 4 * BetheCoefs[9]^3 * (1 - BetheCoefs[9])^1
	Else
	X_C4 = 0
	X_C3 = 0
	X_C2 = 0
	X_C1 = 0
	Endif

	If (BetheCoefs[24] == 8)
	X_B8 = (1 - BetheCoefs[10])^8
	X_B7 = 8 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^7
	X_B6 = 28 * BetheCoefs[10]^2 * (1 - BetheCoefs[10])^6
	X_B5 = 56 * BetheCoefs[10]^3 * (1 - BetheCoefs[10])^5
	X_B4 = 70 * BetheCoefs[10]^4 * (1 - BetheCoefs[10])^4
	X_B3 = 56 * BetheCoefs[10]^5 * (1 - BetheCoefs[10])^3	
	X_B2 = 28 * BetheCoefs[10]^6 * (1 - BetheCoefs[10])^2
	X_B1 = 8 * BetheCoefs[10]^7 * (1 - BetheCoefs[10])^1
	Endif
	
	If (BetheCoefs[24] == 7)
	X_B8 = 0
	X_B7 = 1 * (1 - BetheCoefs[10])^7
	X_B6 = 7 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^6
	X_B5 = 21 * BetheCoefs[10]^2 * (1 - BetheCoefs[10])^5
	X_B4 = 35 * BetheCoefs[10]^3 * (1 - BetheCoefs[10])^4
	X_B3 = 35 * BetheCoefs[10]^4 * (1 - BetheCoefs[10])^3	
	X_B2 = 21 * BetheCoefs[10]^5 * (1 - BetheCoefs[10])^2
	X_B1 = 7 * BetheCoefs[10]^6 * (1 - BetheCoefs[10])^1
	Endif

	If (BetheCoefs[24] == 6)
	X_B8 = 0
	X_B7 = 0
	X_B6 = 1 * (1 - BetheCoefs[10])^6
	X_B5 = 6 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^5
	X_B4 = 15 * BetheCoefs[10]^2 * (1 - BetheCoefs[10])^4
	X_B3 = 20 * BetheCoefs[10]^3 * (1 - BetheCoefs[10])^3	
	X_B2 = 15 * BetheCoefs[10]^4 * (1 - BetheCoefs[10])^2
	X_B1 = 6 * BetheCoefs[10]^5 * (1 - BetheCoefs[10])^1
	Endif

	If (BetheCoefs[24] == 5)
	X_B8 = 0
	X_B7 = 0
	X_B6 = 0
	X_B5 = 1 * (1 - BetheCoefs[10])^5
	X_B4 = 5 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^4
	X_B3 = 10 * BetheCoefs[10]^2 * (1 - BetheCoefs[10])^3	
	X_B2 = 10 * BetheCoefs[10]^3 * (1 - BetheCoefs[10])^2	
	X_B1 = 5 * BetheCoefs[10]^4 * (1 - BetheCoefs[10])^1
	Endif

	If (BetheCoefs[24] == 4)
	X_B8 = 0
	X_B7 = 0
	X_B6 = 0
	X_B5 = 0 
	X_B4 = 1 * (1 - BetheCoefs[10])^4
	X_B3 = 4 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^3	
	X_B2 = 6 * BetheCoefs[10]^2 * (1 - BetheCoefs[10])^2	
	X_B1 = 4 * BetheCoefs[10]^3 * (1 - BetheCoefs[10])^1	
	Endif
	
	If (BetheCoefs[24] == 3)
	X_B8 = 0
	X_B7 = 0
	X_B6 = 0
	X_B5 = 0
	X_B4 = 0
	X_B3 = 1 * (1 - BetheCoefs[10])^3	
	X_B2 = 3 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^2	
	X_B1 = 3 * BetheCoefs[10]^2 * (1 - BetheCoefs[10])^1
	Endif
	
	If (BetheCoefs[24] == 2)
	X_B8 = 0
	X_B7 = 0
	X_B6 = 0
	X_B5 = 0
	X_B4 = 0
	X_B3 = 0
	X_B2 = 1 * (1 - BetheCoefs[10])^2
	X_B1 = 2 * BetheCoefs[10]^1 * (1 - BetheCoefs[10])^1
	Endif
		
	X_A0 = BetheCoefs[9]^(BetheCoefs[23])
	X_B0 = BetheCoefs[10]^(BetheCoefs[24])
	X_C0 = BetheCoefs[9]^(BetheCoefs[23] / 2)
	
	N_A8 = X_A8 * BetheCoefs[20]
	N_A7 = X_A7 * BetheCoefs[20]
	N_A6 = X_A6 * BetheCoefs[20]
	N_A5 = X_A5 * BetheCoefs[20]
	N_A4 = X_A4 * BetheCoefs[20]
	N_A3 = X_A3 * BetheCoefs[20]
	N_A2 = X_A2 * BetheCoefs[20]
	N_A1 = X_A1 * BetheCoefs[20]
	N_A0 = X_A0 * BetheCoefs[20]
	
	N_B8 = X_B8 * BetheCoefs[21]
	N_B7 = X_B7 * BetheCoefs[21]
	N_B6 = X_B6 * BetheCoefs[21]
	N_B5 = X_B5 * BetheCoefs[21]
	N_B4 = X_B4 * BetheCoefs[21]
	N_B3 = X_B3 * BetheCoefs[21]
	N_B2 = X_B2 * BetheCoefs[21]
	N_B1 = X_B1 * BetheCoefs[21]
	N_B0 = X_B0 * BetheCoefs[21]
	
	N_C4 = X_C4 * BetheCoefs[22]
	N_C3 = X_C3 * BetheCoefs[22]
	N_C2 = X_C2 * BetheCoefs[22]
	N_C1 = X_C1 * BetheCoefs[22]
	N_C0 = X_C0 * BetheCoefs[22]
	
	BetheCoefs[28] = N_A8
	BetheCoefs[29] = N_A7
	BetheCoefs[30] = N_A6
	BetheCoefs[31] = N_A5
	BetheCoefs[32] = N_A4
	BetheCoefs[33] = N_A3
	BetheCoefs[34] = N_A2
	BetheCoefs[35] = N_A1
	BetheCoefs[36] = N_A0
	BetheCoefs[37] = N_B8
	BetheCoefs[38] = N_B7
	BetheCoefs[39] = N_B6
	BetheCoefs[40] = N_B5
	BetheCoefs[41] = N_B4
	BetheCoefs[42] = N_B3
	BetheCoefs[43] = N_B2	
	BetheCoefs[44] = N_B1
	BetheCoefs[45] = N_B0	
	BetheCoefs[46] = N_C4	
	BetheCoefs[47] = N_C3
	BetheCoefs[48] = N_C2
	BetheCoefs[49] = N_C1
	BetheCoefs[50] = N_C0
		
	R_A = BetheCoefs[20] / (BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22])
	R_B = BetheCoefs[21] / (BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22])
	R_C = BetheCoefs[22] / (BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22])
	
	BetheCoefs[11] = (X_A8 + X_A7 + X_A6 + X_A5 + X_A4 + X_A3) * R_A + (X_C4 + X_C3) * R_C + (X_B8 + X_B7 + X_B6 + X_B5 + X_B4 + X_B3) * R_B
	BetheCoefs[12] = (4 * X_A8 + 3.5 * X_A7 + 3 * X_A6 + 2.5 * X_A5 + 2 * X_A4 + 1.5 * X_A3) * R_A + (2 * X_C4 + 1.5 * X_C3) * R_C + (4 * X_B8 + 3.5 * X_B7 + 3 * X_B6 + 2.5 * X_B5 + 2 * X_B4 + 1.5 * X_B3) * R_B
	BetheCoefs[13] = BetheCoefs[12] - BetheCoefs[11]
	
//	BetheCoefs[14] = (BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22]) * BetheCoefs[11]
//	BetheCoefs[15] = (BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22]) * BetheCoefs[12]
//	BetheCoefs[16] = (BetheCoefs[20] + BetheCoefs[21] + BetheCoefs[22]) * BetheCoefs[13]
	
	BetheCoefs[14] = (X_A8 + X_A7 + X_A6 + X_A5 + X_A4 + X_A3) * BetheCoefs[20] + (X_B8 + X_B7 + X_B6 + X_B5 + X_B4 + X_B3) * BetheCoefs[21] + (X_C4 + X_C3) * BetheCoefs[22]
	BetheCoefs[15] = (4 * X_A8 + 3.5 * X_A7 + 3 * X_A6 + 2.5 * X_A5 + 2 * X_A4 + 1.5 * X_A3) * BetheCoefs[20] + (4 * X_B8 + 3.5 * X_B7 + 3 * X_B6 + 2.5 * X_B5 + 2 * X_B4 + 1.5 * X_B3) * BetheCoefs[21] + (2 * X_C4 + 1.5 * X_C3) * BetheCoefs[22]
	BetheCoefs[16] = BetheCoefs[15] - BetheCoefs[14]

	BetheCoefs[26] = BetheCoefs[20] * BetheCoefs[17] * X_A0 + BetheCoefs[21] * BetheCoefs[18] * X_B0 + BetheCoefs[22] * BetheCoefs[17] / 2 * X_C0
	BetheCoefs[27] = BetheCoefs[19] - BetheCoefs[26]
	
end


Function BetheCalc()

	wave BetheCoefs
	
	FindRoots/L=0/H=1 BetheFunction, BetheCoefs
	
	return V_root
	
End

Function BetheFunction(w,x)

	Variable x
	Wave w
	
	return (1 - w[5]) + w[5] * ((1 - w[6]) + w[6] * ((1 - w[7]) * x^(w[23] - 1) + w[7] * x^(w[23]/2 - 1)))^(w[24] - 1) - x	
	
End

Function Graph_alpha()

	wave alpha, alpha_xi, alpha_mu, alpha_nu, alpha_cout, BetheCoefs
	wave N_A8, N_A7, N_A6, N_A5, N_A4, N_A3, N_A2, N_A1, N_A0, N_B2, N_B1, N_B0, N_C4, N_C3, N_C2, N_C1, N_C0
	variable i
	
	For (i = 0; i <= 100; i+=1)
	
	BetheCoefs[7] = i / 100
	
	Bethe()
	
	alpha_mu[i] = BetheCoefs[14]
	alpha_nu[i] = BetheCoefs[15]
	alpha_xi[i] = BetheCoefs[16]
	alpha_cout[i] = BetheCoefs[26]
	
	N_A8[i] = BetheCoefs[28]
	N_A7[i] = BetheCoefs[29]
	N_A6[i] = BetheCoefs[30]
	N_A5[i] = BetheCoefs[31]
	N_A4[i] = BetheCoefs[32]
	N_A3[i] = BetheCoefs[33]
	N_A2[i] = BetheCoefs[34]
	N_A1[i] = BetheCoefs[35]
	N_A0[i] = BetheCoefs[36]
	N_B2[i] = BetheCoefs[43]
	N_B1[i] = BetheCoefs[44]
	N_B0[i] = BetheCoefs[45]
	N_C4[i] = BetheCoefs[46]
	N_C3[i] = BetheCoefs[47]
	N_C2[i] = BetheCoefs[48]
	N_C1[i] = BetheCoefs[49]
	N_C0[i] = BetheCoefs[50]	

	If (alpha_mu[i] == 0)
	alpha_mu[i] = -1
	endif
	If (alpha_nu[i] == 0)
	alpha_nu[i] = -1
	endif
	If (alpha_xi[i] == 0)
	alpha_xi[i] = -1
	endif
	If (BetheCoefs[5] > 1)
	alpha_mu[i] = -1
	alpha_nu[i] = -1
	alpha_xi[i] = -1
	endif
	
	endfor

End

Function Graph_connectivity()

	wave probability, probability_xi, probability_mu, probability_nu, probability_cout, BetheCoefs
	variable i
	
	For (i = 0; i <= 100; i+=1)
	
	BetheCoefs[4] = i / 100
	
	Bethe()
	
	probability_mu[i] = BetheCoefs[14]
	probability_nu[i] = BetheCoefs[15]
	probability_xi[i] = BetheCoefs[16]
	probability_cout[i] = BetheCoefs[26]
	
	If (probability_mu[i] == 0)
	probability_mu[i] = -1
	endif
	If (probability_nu[i] == 0)
	probability_nu[i] = -1
	endif
	If (probability_xi[i] == 0)
	probability_xi[i] = -1
	endif
	If (BetheCoefs[5] > 1)
	probability_mu[i] = -1
	probability_nu[i] = -1
	probability_xi[i] = -1
	endif
		
	endfor

End

Function Calculation_alpha()
	
	variable p, xi, alpha
	variable i, j
	wave BetheCoefs, calc_xi, Input_xi, Output_alpha
	
	p = 1

For (j = 0; j <= 6; j+=1)

	xi = Input_xi[j]

For (i = 0; i <= 1000; i+=1)
	
	BetheCoefs[3] = p
	BetheCoefs[7] = i / 1000
	Bethe()
	calc_xi[i] = (BetheCoefs[16] - xi)^2
	
	If (calc_xi[i] == 0)
	calc_xi[i] = -1
	endif
	If (BetheCoefs[5] > 1)
	calc_xi[i] = -1
	endif
		
endfor

	wavestats calc_xi
	Output_alpha[j] = V_minloc / 1000

endfor
	
End

Function Calculation_pc()

	wave alpha_pc, BetheCoefs
	variable i
	
	p_net()
	p_A()
	p_B()
	
	For (i = 0; i <= 100; i+=1)
	
	BetheCoefs[7] = i / 100
	
	alpha_pc[i] = (4 * BetheCoefs[4] * (1 - BetheCoefs[4]) / ((BetheCoefs[24] - 1) * (BetheCoefs[23] - 1 - BetheCoefs[23] / 2 * BetheCoefs[7])))^0.5
	
	Endfor
	
End

Function Calculation_rc()

	wave alpha_rc, BetheCoefs
	variable i
	
	p_net()
	p_A()
	p_B()
	
	For (i = 0; i <= 100; i+=1)
	
	BetheCoefs[7] = i / 100
	
	alpha_rc[i] = (BetheCoefs[24] - 1) * (BetheCoefs[23] - 1 - BetheCoefs[23] / 2 * BetheCoefs[7]) / (BetheCoefs[2]^2 + (BetheCoefs[24] - 1) * (BetheCoefs[23] - 1 - BetheCoefs[23] / 2 * BetheCoefs[7]))
	
	Endfor
	
End