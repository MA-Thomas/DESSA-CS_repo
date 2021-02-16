package main

import (
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"
)

type Accuracy float64

func Prod_ExpErfc(a, b float64) *big.Float {
	//fmt.Println("ExpErfc_TaylorApprox file. a,b are: ", a, b)
	Product := new(big.Float).SetFloat64(0.0)

	if math.IsInf(math.Erfc(b), 0) {
		panic("Erfc should never be infinite")
	}


	if math.Erfc(b) == 0 || b > 27.0 {
		// This case is to handle when erfc(b) underflows to 0.
		// This occurs around 27.22353386859086 < b < 27.22647234390665.
		// However, the accuracy of math.Erfc(b) seems to go down somewhere around
		// b = 27.08801942429744. Maybe earlier. Therefore, this code applies as
		// soon as b crosses the 27.000... threshold.

		/* Product = exp(a)erfc(b)
							 = exp(b)erfc(b)*exp(a-b)

		We can make use of the Taylor Expansion at infinity for  exp(b)erfc(b)*exp(a-b)
		Note: in wolframalpha, enter: Series[E^b Erfc[b] E^[a-b] , {b, Infinity, 20}]
		https://www.wolframalpha.com/input/?i=Series%5BE%5Eb+Erfc%5Bb%5D+E%5E%5Ba-b%5D+%2C+%7Bb%2C+Infinity%2C+20%7D%5D



		Now take the log.



		At the end, exponentiate this.

		*/
		//A := new(big.Float).SetFloat64(a)
		B := new(big.Float).SetFloat64(b)
		sqrtPi := new(big.Float).SetFloat64( math.Sqrt(math.Pi) )

		order1 := new(big.Float).SetFloat64( 1.0 )
		order1.Quo(order1, sqrtPi)
		order1.Quo(order1,B)


		// log( 1/[2*sqrt(pi)*b^3] )
		order3 := new(big.Float).SetFloat64( 1.0 )
		order3.Quo(order3, sqrtPi)
		order3.Quo(order3,new(big.Float).SetFloat64(2.0))
		order3.Quo(order3,B)
		order3.Quo(order3,B)
		order3.Quo(order3,B)


		// log( 3/[4*sqrt(pi)*b^5] )
		order5 := new(big.Float).SetFloat64( 3.0 )
		order5.Quo(order5, sqrtPi)
		order5.Quo(order5,new(big.Float).SetFloat64(4.0))
		order5.Quo(order5,B)
		order5.Quo(order5,B)
		order5.Quo(order5,B)
		order5.Quo(order5,B)
		order5.Quo(order5,B)


		// log( 15/[8*sqrt(pi)*b^7] )
		order7 := new(big.Float).SetFloat64( 15.0 )
		order7.Quo(order7, sqrtPi)
		order7.Quo(order7,new(big.Float).SetFloat64(8.0))
		order7.Quo(order7,B)
		order7.Quo(order7,B)
		order7.Quo(order7,B)
		order7.Quo(order7,B)
		order7.Quo(order7,B)
		order7.Quo(order7,B)
		order7.Quo(order7,B)


		// log( 105/[16*sqrt(pi)*b^9] )
		order9 := new(big.Float).SetFloat64( 105.0 )
		order9.Quo(order9, sqrtPi)
		order9.Quo(order9,new(big.Float).SetFloat64(16.0))
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)
		order9.Quo(order9,B)


		// log( 945/[32*sqrt(pi)*b^11] )
		order11 := new(big.Float).SetFloat64( 945.0 )
		order11.Quo(order11, sqrtPi)
		order11.Quo(order11,new(big.Float).SetFloat64(32.0))
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)
		order11.Quo(order11,B)


		// log( 10395/[64*sqrt(pi)*b^13] )
		order13 := new(big.Float).SetFloat64( 10395.0 )
		order13.Quo(order13, sqrtPi)
		order13.Quo(order13,new(big.Float).SetFloat64(64.0))
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)
		order13.Quo(order13,B)


		// log( 135135/[128*sqrt(pi)*b^15] )
		order15 := new(big.Float).SetFloat64( 135135.0 )
		order15.Quo(order15, sqrtPi)
		order15.Quo(order15,new(big.Float).SetFloat64(128.0))
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)
		order15.Quo(order15,B)

		order17 := new(big.Float).SetFloat64( 2027025.0 )
		order17.Quo(order17, sqrtPi)
		order17.Quo(order17,new(big.Float).SetFloat64(256.0))
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)
		order17.Quo(order17,B)

		Sum := new(big.Float).SetFloat64(0.0)
		Sum.Add(Sum,order1)
		Sum.Sub(Sum,order3)
		Sum.Add(Sum,order5)
		Sum.Sub(Sum,order7)
		Sum.Add(Sum,order9)
		Sum.Sub(Sum,order11)
		Sum.Add(Sum,order13)
		Sum.Sub(Sum,order15)
		Sum.Add(Sum,order17)
		
		SumFloat64, _ := Sum.Float64()
		logSum := a-math.Pow(b,2) + math.Log(SumFloat64)

		Product = new(big.Float).SetFloat64(math.Exp(logSum))
		//fmt.Println("a,b,logSum,product: ",a,b,logSum,Product)


		//fmt.Println("ExpErfc_TaylorApprox file, erfc(b)==0 CASE. b = ")
		//Product = new(big.Float).SetFloat64(math.Exp(-math.Pow(b, 2.0) + math.Log((4*math.Pow(b, 4.0)-2*math.Pow(b, 2.0)+3)/(4*math.Sqrt(math.Pi)*math.Pow(b, 5.0))) + a))
		return Product
	}

	if math.IsInf(math.Exp(a), 0) {
		// -------------------------------------------------------------------------
		// -------------------------------------------------------------------------
		// NOTE: this method seems great, but numerical results don't agree with
		// WolframAlpha.
		/*
			    // Main idea: exp(a)erfc(b) => exp(b)erfc(b) * exp(a-b)
					// This makes the computation easier because exp(b)erfc(b) will be a very
					// small number and exp(a-b) is less likely to be +inf as compared with
					// exp(a).
					// Eventually, even exp(a-b) will become +inf and again the code will break.
					// To put this off, we can simply replace exp(a-b) with a set of equivalent
					// multiplications: Prod over h: exp[(a-b)/h]
					Exp1 := new(big.Float).SetFloat64(math.Exp(b))
					Erfc := new(big.Float).SetFloat64(math.Erfc(b))
					//Erfc := new(big.Float).SetFloat64(1.0 - math.Erf(b))
					Product.Mul(Exp1, Erfc)
					fmt.Println("a,b and exp(b)erfc(b) are ", a, b, Product)
					fmt.Println(a - b)
					remainderExponent := a - b
					subdivision := 50 // exp(a-b) = { exp[(a-b)/subdivision)]*exp[(a-b)/subdivision)]*...*exp[(a-b)/subdivision)] }
					for h := 0; h < subdivision; h++ {
						exponent := float64(remainderExponent) / float64(subdivision)
						Product.Mul(Product, new(big.Float).SetFloat64(math.Exp(exponent)))
						//fmt.Println("Full product is ", Product)
					}

					//fmt.Println("remainder Exp is ", math.Exp(a-b))
					//fmt.Println("Full product is ", Product)
		*/

		/*
		   Something strange. For the following inputs A,B, the product of exp(A)erfc(B) is:
		   A,B,Product:  738.7207963063141 27.179418616046853 0.01969252521720255

		   According to wolfram alpha, the product should be
		   0.020743949493643551070503970098523...

		*/
		// -------------------------------------------------------------------------
		// -------------------------------------------------------------------------

		// -------------------------------------------------------------------------
		// -------------------------------------------------------------------------
		// ----------- TRY EVALUATING log(exp(a)*erfc(b)), then exponentiating
		/*
					    log(exp(a)*erfc(b)) = a + log(erfc(b))
					    NOTE: erfc(b) will, in scientific notation, look like 4.3324...e-323
					    Golang will fail to evaluate something like log(erfc(4.3324...e-323)) accurately.
					    The discrepancy is not negligable.
					    Solution: convert log(erfc(4.3324...e-323)) => log(4.3324...) - log(10)log10(e323)
					                                                = log(4.3324...) - log(10)*323

					   FINAL RESULT: exp(a)*erfc(b) = exp( a + log(significant digits) - log(10)*exponent )

						 Note: Normally to get the exponent of N, one might try floor(log10(N)),
			       however, this leads to the same issue of Golang's poor accuracy.
			       I had to resort to using string maniplulations instead.

						 Note: This formulation assumes erfc(b) does not underflow to 0.
						       Through tests, I discovered that the accuracy of math.Erfc(b)
									 - for values of y leading up to the underflow value of y -
									 declines a bit.
									 Therefore, we don't want to go to the Taylor approximation
									 case earlier in this file exactly when math.Erfc(b)=0. We
									 want to go to that approximation at slightly smaller b values.


		*/

		//fmt.Println("NAIVE: a,b, erfc(b), log(erfc(b)), LOGofProd, Product are: ", a, b, math.Erfc(b), math.Log(math.Erfc(b)), a+math.Log(10.0)*math.Log10(math.Erfc(b)), math.Exp(a+math.Log(10.0)*math.Log10(math.Erfc(b))))

		ERFCb := math.Erfc(b)
		S := fmt.Sprintf("%.16e", ERFCb)
		//fmt.Println("S is: ", S)
		//runeS := []rune(S)
		exponent := ""
		substring := ""
		significantDigits := ""
		for i := 0; i < len(S); i++ {
			substring = S[i : i+1]
			if strings.Contains(substring, "e") {
				exponent = S[i+2 : len(S)]
				significantDigits = S[:i]
			}
		}
		SigDigits, _ := strconv.ParseFloat(significantDigits, 64)
		Exponent, _ := strconv.ParseFloat(exponent, 64)

		//fmt.Println("significantDigits and exponent are: ", SigDigits, Exponent)
		//fmt.Println("log(significantDigits), log(exponent), difference are: ", math.Log(SigDigits), math.Log(10)*Exponent, math.Log(SigDigits)-math.Log(10)*Exponent)

		LOGofProd := a + math.Log(SigDigits) - math.Log(10)*Exponent
		//fmt.Println("TRUE: Log( exp(x)*erfc(y) ) is: ", LOGofProd)
		ExpLOGofProd := math.Exp(LOGofProd)
		//fmt.Println("TRUE: ExpLOGofProd (i.e. Product) is: ", ExpLOGofProd)
		Product = new(big.Float).SetFloat64(ExpLOGofProd)

	} else {

		Exp := new(big.Float).SetFloat64(math.Exp(a))
		Erfc := new(big.Float).SetFloat64(math.Erfc(b))
		Product.Mul(Exp, Erfc)

	}

	return Product
}

// MT - from https://steemit.com/tutorial/@gopher23/power-and-root-functions-using-big-float-in-golang
func Pow(a *big.Float, e int64) *big.Float {
	result := new(big.Float).SetFloat64(1.0)

	if e == int64(0) {
		// return 1

	} else if e == int64(1) {
		result.Copy(a)

	} else if e == int64(-1) {
		result.Copy(a)

		result.Quo(new(big.Float).SetFloat64(1.0), result)

	} else {
		result.Copy(a)
		for i := int64(1); i < e; i++ {
			result.Mul(result, a)
		}
	}
	return result
}

// MT - from https://play.golang.org/p/53TmmygltkR
func factorial(x *big.Int) *big.Int {
	n := big.NewInt(1)
	if x.Cmp(big.NewInt(0)) == 0 {
		return n
	}
	return n.Mul(x, factorial(n.Sub(x, n)))
}
