package main

import (

  "math"
  "math/big"
  _"fmt"

)
type Accuracy float64

func Prod_ExpErfc(a, b float64) *big.Float {
  //var val float64 = 0.0
  Product := new(big.Float).SetFloat64(0.0)

  if math.IsInf(math.Erfc(b), 0) {
    panic("Erfc should never be infinite")
  }

  if math.Erfc(b) == 0 {
    return Product
  }

  if math.IsInf(math.Exp(a), 0) {

      //fmt.Println("exp term IS infinite")
      Erfc := new(big.Float).SetFloat64(math.Erfc(b))

      A := new(big.Float).SetFloat64(a)
      //B := new(big.Float).SetFloat64(b)

      One := new(big.Float).SetFloat64(1.0) // used to turn erf into erfc
      //MinusOne := new(big.Float).SetFloat64(-1.0)

      // Taylor approx of exp(x) starts with 1 + x (0-th and 1st terms)
      Taylor_Approx_Exp := new(big.Float).SetFloat64(a)
      Taylor_Approx_Exp.Add(Taylor_Approx_Exp,One)

      // Taylor approx of erf(z) starts with z - z^3/3 (0-th and 1-st terms)
      Taylor_Approx_Erf := new(big.Float).SetFloat64(b - math.Pow(b,3)/3.0)

      numTerms := 100

      factorialTerm := new(big.Float).SetFloat64(1.0)
      nextExpNumerator := new(big.Float).SetFloat64(a)

      Temp_Exp := new(big.Float).SetFloat64(0.0)

      Product.Mul(Taylor_Approx_Exp, Taylor_Approx_Erf)

      Exp_converged := false
      Exp_convergence_delta := 1e-7

      // i is the index in the Taylor series sum
      for i := 2; i <= numTerms; i++ {
          //fmt.Println("Taylor approx of erf is ", Taylor_Approx_Erf)
          I := new(big.Float).SetFloat64(float64(i))
          factorialTerm.Mul(factorialTerm, I) // implementing i-th factorial

          // Exp
          if Exp_converged == false {
            nextExpNumerator.Mul(nextExpNumerator,A) // implementing raise to i-th power
            Temp_Exp.Quo(nextExpNumerator,factorialTerm)

            Exp_ratio := new(big.Float).SetFloat64(1.0)
            Exp_ratio.Mul(Taylor_Approx_Exp, One)
            Taylor_Approx_Exp.Add(Taylor_Approx_Exp, Temp_Exp) // update exp approximation
            Exp_ratio.Quo(Taylor_Approx_Exp, Exp_ratio)// update the ratio of current approx to prev approx
            //fmt.Println("Taylor approx of exp: ",Taylor_Approx_Exp)
            //fmt.Println("Exp_ratio is ",Exp_ratio)

            check,_ := Exp_ratio.Float64()
            //fmt.Println("1 - exp ratio is ",math.Abs(1-check))
            if math.Abs(1-check) < Exp_convergence_delta {
              //fmt.Println("exp converged at i = ",i)
              Exp_converged = true
            }
          }

      }

      Product.Mul(Taylor_Approx_Exp,Erfc)

      //fmt.Println("Taylor approx of exp is ", Taylor_Approx_Exp)
      //fmt.Println("Erfc is ",Erfc)
      //fmt.Println("Product is ",Product)





  } else {



    //fmt.Println("exp term is not infinite")
    Exp := new(big.Float).SetFloat64(math.Exp(a))
    Erfc := new(big.Float).SetFloat64(math.Erfc(b))
    Product.Mul(Exp, Erfc)

    //fmt.Println("Exp is ",Exp)
    //fmt.Println("Erfc is ",Erfc)

  }

  return Product
}
