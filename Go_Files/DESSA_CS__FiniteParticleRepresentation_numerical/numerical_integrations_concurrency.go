package main

import (
  "fmt"
  "time"
  "math"
  _"math/rand"
  "math/big"
  "os"
  "log"
  "encoding/csv"

)

func main() {
/*

  // TESTING  Prod_ExpErfc()
  startTime := time.Now()
  g := Prod_ExpErfc(-1e-6, 5)
  elapsedTime := time.Since(startTime)
  fmt.Println("product = ",g)
  fmt.Println("run time = ",elapsedTime)
*/
///*

  startTime := time.Now()


  var contactR float64
  contactR = 0.01

  var n_sigma float64
  n_sigma = 5.0

  var Dmonomer float64
  Dmonomer = 1.0

  var T_hardlimit_diffuse float64
  T_hardlimit_diffuse = 40.0

  var earliestPropensityTime float64
  earliestPropensityTime = 1.0e-6

  // the intrinsic reaction rate is used to compute the PDF (in the finite particle size representation)
  // It is used *again* to evaluate the integrated reaction propensity from the
  // integrated PDF
  var c float64
  c = 5.0e-1

  numDistances := 3000
  dmin := contactR
  dmax := 2*n_sigma*math.Sqrt(6*Dmonomer*T_hardlimit_diffuse)
  distanceSpacing := (dmax - dmin) / float64(numDistances)

  dist_list := make([]float64,numDistances)
  dist_list[0] = dmin
  for i := 1; i < numDistances; i++ {
    dist_list[i] = dist_list[i-1] + distanceSpacing
  }

  numTimes1 := 40000 // i.e. # time points within [0,1]
  numTimes2 := 10000 // i.e. # time points within [1,T_hardlimit_diffuse]
  tmin := earliestPropensityTime
  tmax := T_hardlimit_diffuse

  shortTime := 1e-2
  timeSpacing1 := (shortTime - tmin) / float64(numTimes1)
  timeSpacing2 := (tmax - shortTime) / float64(numTimes2)

  times_list := make([]float64,numTimes1 + numTimes2)
  times_list[0] = tmin
  for i := 1; i < (numTimes1 + numTimes2); i++ {
    if times_list[i-1] + timeSpacing1 < shortTime {
    times_list[i] = times_list[i-1] + timeSpacing1

    }else {
      times_list[i] = times_list[i-1] + timeSpacing2
    }
  }


  //fmt.Println("times_list is: \n", times_list)
  //fmt.Println("timeSpacing1 is: \n", timeSpacing1)
  //fmt.Println("timeSpacing2 is: \n", timeSpacing2)
  //fmt.Println("var_list is: \n",var_list)
  //fmt.Println("lengths dist_list,time_list: ",len(dist_list),len(times_list))

  PDFs, Integrated_PDFs := IntegrateN(contactR, 2*Dmonomer, c, dist_list, times_list)


  // The integrated propensities are just the integrated PDFs multiplied by the
  // intrinsic rate constant
  Integrated_Propensities := make([][]float64,len(Integrated_PDFs))
  for i := 0; i < len(Integrated_PDFs); i++ {
    Integrated_Propensities[i] = make([]float64, len(Integrated_PDFs[i]))
    for j := 0; j < len(Integrated_Propensities[i]); j++ {
      Integrated_Propensities[i][j] = c * Integrated_PDFs[i][j]
    }
  }

  // Run time to perform integrations
  elapsedTime := time.Since(startTime)
  fmt.Println("run time (integrations): ",elapsedTime)

  // Compute trapezoidal error bounds.
  // Start by computing the max second derivative for each distance value
  startTime = time.Now()
  maxSecondDerivatives := make([]float64,numDistances)

  for i, Row := range PDFs {
    maxSecondDerivatives[i] = maxSecondDerivative(Row, timeSpacing1, timeSpacing2, numTimes1, numTimes2)
  }

  errorBounds := make([][]float64,len(dist_list) )
  for i := 0; i < len(errorBounds); i++ {
    errorBounds[i] = make([]float64, len(times_list))

    for N, timeInSec := range times_list {

      // This N appears in the denom of the error formula. We cannot have denom==0
      if N == 0 || N == numTimes1 {
        continue
      }

      if N < numTimes1 { // error bounds for N = [1,numTimes1)
        errorBounds[i][N] = math.Pow((timeInSec - tmin),3) * maxSecondDerivatives[i] / (12*math.Pow(float64(N),2))
      } else { // error bounds for N = [numTimes1+1:)
        t_start := times_list[numTimes1-1]
        errorBounds[i][N] = errorBounds[i][numTimes1-1] + math.Pow((timeInSec - t_start),3) * maxSecondDerivatives[i] / (12*math.Pow(float64(N-numTimes1),2))
      }
    }
  }

  // Run time to eval max second derivatives at each distance
  elapsedTime = time.Since(startTime)
  fmt.Println("run time (second derivatives/error bounds): ",elapsedTime)
  fmt.Println("There are ", numTimes1, " time points between 0 and ", shortTime)
  fmt.Println("There are ", numTimes2, " time points between ", shortTime, "and ", tmax)
  fmt.Println("All done with computations. Next, save to files...")
  //fmt.Println("max 2nd derivatives are: \n",maxSecondDerivatives)
  //fmt.Println("Error bounds d_ind=0 are: \n",errorBounds[0])



//*/

  // Write outputs to files
  fmt.Println("Writing output to files")
// /*
  fileName := "integratedPropensities_FINITE.csv"
  success := writeCSV(Integrated_Propensities, fileName)
  fmt.Println(success)

  fileName := "integratedPDFs_FINITE.csv"
  success := writeCSV(Integrated_PDFs, fileName)
  fmt.Println(success)

  fileName = "PDFs_FINITE.csv"
  success = writeCSV(PDFs, fileName)
  fmt.Println(success)

  fileName = "maxSecondDerivatives_FINITE.csv"
  success = writeCSV2(maxSecondDerivatives, fileName)
  fmt.Println(success)

  fileName = "errorBounds_FINITE.csv"
  success = writeCSV(errorBounds, fileName)
  fmt.Println(success)

  fileName = "distance_list_FINITE.csv"
  success = writeCSV2(dist_list, fileName)
  fmt.Println(success)

  fileName = "times_list_FINITE.csv"
  success = writeCSV2(times_list, fileName)
  fmt.Println(success)

// */


  return
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

type integrationOutput struct {
  PDFs, integratedPDFs []float64
  d_index int
}

func IntegrateN(contactR, D, c float64, dist_list, times_list []float64) ([][]float64,[][]float64) {

  // Initialize the PDF and Integrated PDF matrices
  PDF_arrays := make([][]float64, len(dist_list))
  for i := 0; i < len(PDF_arrays); i++ {
    PDF_arrays[i] = make([]float64, len(times_list))
  }
  Integrated_PDF_arrays := make([][]float64, len(dist_list))
  for i := 0; i < len(Integrated_PDF_arrays); i++ {
    Integrated_PDF_arrays[i] = make([]float64, len(times_list))
  }

  numDistances := len(dist_list)
  numTimes := len(times_list)



  channel := make(chan integrationOutput,numDistances)
  for d_ind := 0; d_ind < numDistances; d_ind++ {
  //for d_ind := 0; d_ind < 1; d_ind++ {
    //nan_reached = 0
    Distance := dist_list[d_ind]
    if d_ind%100 == 0 {
      fmt.Println("Distance, d_ind = ",Distance,d_ind)
    }

    //channel := make(chan integrationOutput)
    go evalIntegralN(Distance,contactR,D,c,times_list,numTimes,d_ind,channel)
    /*
    // With this code structure, output is not set until the channel is ready to
    // send. Hence we get no performance gain bc the code is still running
    // entirely sequentially. Instead, we should only keep the go routine calls
    // in this for loop, and then in a subsequent for loop, receive the results
    // as they finish.
    output := <- channel
    index := output.d_index
    PDF_arrays[index] = output.PDFs
    Integrated_PDF_arrays[index] = output.integratedPDFs
    //fmt.Println("Index is :",index)
    //fmt.Println("channel output is :\n", output)
    */


  }


  for i := 1; i <= numDistances; i++ {
    output := <- channel // this line will execute as soon as the channel has something to send
    index := output.d_index
    PDF_arrays[index] = output.PDFs
    Integrated_PDF_arrays[index] = output.integratedPDFs

    //fmt.Println("channel output is :\n", output)
    fmt.Println("d_ind ",index, " completed.")
  }
  close(channel)
  return PDF_arrays, Integrated_PDF_arrays
}



func evalIntegralN(r0,contactR,D,c float64, times_list []float64, numTimes, d_ind int, Ch chan integrationOutput) {
  var Output integrationOutput
  probabilityDensity_vector := make([]float64,numTimes)
  Integrated_probabilityDensity_vector := make([]float64,numTimes)

  var B float64 = (1.0 + (c / (4*math.Pi*contactR*D))) / contactR

  result := new(big.Float).SetFloat64(1.0)

  for t_ind := 0; t_ind < numTimes; t_ind++{
    tau := times_list[t_ind]

    term1 := new(big.Float).SetFloat64( (1.0 / (8*math.Pi*contactR*r0)) * (1.0 / math.Sqrt(D*tau*math.Pi)) )
    term2 := new(big.Float).SetFloat64( math.Exp(-math.Pow(contactR-r0, 2)/(4.0*D*tau) ) )
    term3 := new(big.Float).SetFloat64( math.Exp(-math.Pow(r0-contactR, 2)/(4.0*D*tau) ) )

    term4 := new(big.Float).SetFloat64(1.0)
    term4_pre2 := Prod_ExpErfc(math.Pow(B,2)*D*tau + B*(r0-contactR),
    (r0-contactR)/(2*math.Sqrt(D*tau)) + B*math.Sqrt(D*tau))
    term4_pre1 := new(big.Float).SetFloat64( -(2.0*B*math.Sqrt(D*tau*math.Pi)) )

    term4.Mul(term4_pre1, term4_pre2)

    summand := new(big.Float).SetFloat64(0.0)
    summand.Add(term2, term3)
    summand.Add(summand, term4)

    result.Mul(term1, summand)
    probabilityDensity_vector[t_ind], _ = result.Float64()




///*
    // Trapezoidal Integration
    if t_ind > 0 {
      //fmt.Println(times_list[v_ind],times_list[v_ind-1])
      //fmt.Println("integral is ",(times_list[v_ind]-times_list[v_ind-1])*(PDF_arrays[d_ind][v_ind] + PDF_arrays[d_ind][v_ind-1])/2)
      Integrated_probabilityDensity_vector[t_ind] = Integrated_probabilityDensity_vector[t_ind-1] + (times_list[t_ind]-times_list[t_ind-1])*(probabilityDensity_vector[t_ind] + probabilityDensity_vector[t_ind-1])/2
      //fmt.Println("Integrated_probabilityDensity_vector[t_ind]",Integrated_probabilityDensity_vector[t_ind])
    }
//*/
  }

  // Gather the results into a struct
  // Send the struct over the channel, Ch
  // Don't close the chnnel here.
  Output.PDFs = probabilityDensity_vector
  Output.integratedPDFs = Integrated_probabilityDensity_vector
  Output.d_index = d_ind
  Ch <- Output

}



func writeCSV( data [][]float64, fileName string) bool {
  file, err := os.Create(fileName)
  checkError("Cannot create file", err)
  defer file.Close()

  writer := csv.NewWriter(file)
  defer writer.Flush()

  for row, _ := range data {
      str := []string{}
      for _, element := range data[row] {
        text := fmt.Sprintf("%g",element)
        str = append(str, text)
      }

      err := writer.Write(str)
      checkError("Cannot write to file", err)
  }

  return true
}
func checkError(message string, err error) {
    if err != nil {
        log.Fatal(message, err)
    }
}

func writeCSV2( data []float64, fileName string) bool {
  file, err := os.Create(fileName)
  checkError("Cannot create file", err)
  defer file.Close()

  writer := csv.NewWriter(file)
  defer writer.Flush()


  str := []string{}
  for _, element := range data {
    text := fmt.Sprintf("%g",element)
    str = append(str, text)
  }

  err = writer.Write(str)
  checkError("Cannot write to file", err)


  return true
}


func MAXelem(f []float64) (int, float64) {
  max := 0.0
  index := 0
  for i,val := range f {
    if val > max {
      max = val
      index = i
    }
  }
  return index, max
}

func maxSecondDerivative(f []float64, h1, h2 float64, numTimePoints1, numTimePoints2 int) float64 {
  if numTimePoints1 + numTimePoints2 != len(f) {
    panic("len(numTimePoints1) + len(numTimePoints2) should equal len(f)")
  }

  // The time point spacing used to evaluate f (i.e. the PDFs) is constant for
  // the first <numTimePoints1> points, then changes to a new constant value for
  // the remaining <numTimePoints2> points.
  // These spacings are, respectively, h1 and h2.
  // Thus, when computing numerical second derivatives, we need to use the
  // appropriate values.

  f1 := f[:numTimePoints1]
  f2 := f[numTimePoints1:]
  if len(f1) != numTimePoints1 || len(f2) != numTimePoints2 {
    panic("we need: len(f1) == numTimePoints1 && len(f2) == numTimePoints2")
  }


  DDf1 := make([]float64, len(f1))
  for n := 1; n < len(f1)-1; n++ {
    if f1[n-1] == 0 {
      continue
    } else {
      DDf1[n] = ( f1[n+1] - 2.0*f1[n] + f1[n-1] ) / math.Pow(h1,2)
    }
  }

  DDf2 := make([]float64, len(f2))
  for n := 1; n < len(f2)-1; n++ {
    if f2[n-1] == 0 {
      continue
    } else {
      DDf2[n] = ( f2[n+1] - 2.0*f2[n] + f2[n-1] ) / math.Pow(h2,2)
    }
  }

  DDf := append(DDf1,DDf2...)
  _,maxDDf := MAXelem(DDf)
  return maxDDf
}
