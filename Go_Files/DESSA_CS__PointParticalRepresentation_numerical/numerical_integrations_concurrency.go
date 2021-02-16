package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"

	"encoding/csv"
	"log"
	"os"
)

func main() {
	startTime := time.Now()

	var contactR float64
	contactR = 0.01

	var n_sigma float64
	n_sigma = 5.0

	var Dmonomer float64
	Dmonomer = 1.0

	var T_hardlimit_diffuse float64
	T_hardlimit_diffuse = 40 //32 //16.0
	//T_hardlimit_diffuse = 2 // for rebinding time tests

	var earliestPropensityTime float64
	earliestPropensityTime = 1.0e-6

	numDistances := 3000
	dmin := contactR
	dmax := 2 * n_sigma * math.Sqrt(6*Dmonomer*T_hardlimit_diffuse)
	distanceSpacing := (dmax - dmin) / float64(numDistances)

	dist_list := make([]float64, numDistances)
	dist_list[0] = dmin
	for i := 1; i < numDistances; i++ {
		dist_list[i] = dist_list[i-1] + distanceSpacing
	}

	numVariances1 := 40000 // i.e. # time points within [0,1]
	numVariances2 := 10000 // i.e. # time points within [1,T_hardlimit_diffuse]
	//numVariances1 := 100000 // for rebinding time tests
	//numVariances2 := 20000  // for rebinding time tests
	tmin := earliestPropensityTime
	tmax := T_hardlimit_diffuse

	shortTime := 1e-2
	//shortTime := 1e-1 // for rebinding time tests
	timeSpacing1 := (shortTime - tmin) / float64(numVariances1)
	timeSpacing2 := (tmax - shortTime) / float64(numVariances2)

	times_list := make([]float64, numVariances1+numVariances2)
	times_list[0] = tmin
	for i := 1; i < (numVariances1 + numVariances2); i++ {
		if times_list[i-1]+timeSpacing1 < shortTime {
			times_list[i] = times_list[i-1] + timeSpacing1

		} else {
			times_list[i] = times_list[i-1] + timeSpacing2
		}
	}

	var_list := make([]float64, numVariances1+numVariances2)
	Da := Dmonomer
	Db := Dmonomer
	for i := 0; i < len(var_list); i++ {
		var_list[i] = 6*Da*times_list[i] + 6*Db*times_list[i]
	}

	CDFs, Integrated_CDFs := IntegrateN(contactR, dist_list, var_list, times_list)

	fmt.Println("Length of Integrated_CDFs = ", len(Integrated_CDFs))
	fmt.Println("Length CDFs is ", len(CDFs))

	// Run time to perform integrations
	elapsedTime := time.Since(startTime)
	fmt.Println("run time (integrations): ", elapsedTime)

	// Compute trapezoidal error bounds.
	// Start by computing the max second derivative for each distance value
	startTime = time.Now()
	maxSecondDerivatives := make([]float64, numDistances)
	for i, Row := range CDFs {
		maxSecondDerivatives[i] = maxSecondDerivative(Row, timeSpacing1, timeSpacing2, numVariances1, numVariances2)
	}

	errorBounds := make([][]float64, len(dist_list))
	for i := 0; i < len(errorBounds); i++ {
		errorBounds[i] = make([]float64, len(times_list))

		for N, timeInSec := range times_list {

			// This N appears in the denom of the error formula. We cannot have denom==0
			if N == 0 || N == numVariances1 {
				continue
			}

			if N < numVariances1 { // error bounds for N = [1,numVariances1)
				errorBounds[i][N] = math.Pow((timeInSec-tmin), 3) * maxSecondDerivatives[i] / (12 * math.Pow(float64(N), 2))
			} else { // error bounds for N = [numVariances1+1:)
				t_start := times_list[numVariances1-1]
				errorBounds[i][N] = errorBounds[i][numVariances1-1] + math.Pow((timeInSec-t_start), 3)*maxSecondDerivatives[i]/(12*math.Pow(float64(N-numVariances1), 2))
			}
		}
	}

	// Run time to eval max second derivatives at each distance
	elapsedTime = time.Since(startTime)
	fmt.Println("run time (second derivatives/error bounds): ", elapsedTime)
	fmt.Println("There are ", numVariances1, " time points (& variances) between 0 and ", shortTime)
	fmt.Println("There are ", numVariances2, " time points (& variances) between ", shortTime, "and ", tmax)
	fmt.Println("All done with computations. Next, save to files...")

	// /*
	// Write outputs to files
	fmt.Println("Writing output to files")
	fileName := "integratedCDFs.csv"
	success := writeCSV(Integrated_CDFs, fileName)
	fmt.Println(success)

	fileName = "CDFs.csv"
	success = writeCSV(CDFs, fileName)
	fmt.Println(success)

	//fileName = "maxSecondDerivatives.csv"
	//success = writeCSV2(maxSecondDerivatives, fileName)
	//fmt.Println(success)

	//fileName = "errorBounds.csv"
	//success = writeCSV(errorBounds, fileName)
	//fmt.Println(success)

	fileName = "distance_list.csv"
	success = writeCSV2(dist_list, fileName)
	fmt.Println(success)

	fileName = "times_list.csv"
	success = writeCSV2(times_list, fileName)
	fmt.Println(success)

	fileName = "variance_list.csv"
	success = writeCSV2(var_list, fileName)
	fmt.Println(success)
	// */

	return
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

type integrationOutput struct {
	CDFs, integratedCDFs []float64
	d_index              int
}

func IntegrateN(contactR float64, dist_list, var_list, times_list []float64) ([][]float64, [][]float64) {

	// Initialize the CDF and Integrated CDF matrices
	CDF_arrays := make([][]float64, len(dist_list))
	for i := 0; i < len(CDF_arrays); i++ {
		CDF_arrays[i] = make([]float64, len(times_list))
	}
	Integrated_CDF_arrays := make([][]float64, len(dist_list))
	for i := 0; i < len(Integrated_CDF_arrays); i++ {
		Integrated_CDF_arrays[i] = make([]float64, len(times_list))
	}

	numDistances := len(dist_list)
	numVariances := len(var_list)

	channel := make(chan integrationOutput, numDistances)
	for d_ind := 0; d_ind < numDistances; d_ind++ {
		//for d_ind := 0; d_ind < 1; d_ind++ {
		//nan_reached = 0
		Distance := dist_list[d_ind]
		if d_ind%100 == 0 {
			fmt.Println("Distance, d_ind = ", Distance, d_ind)
		}

		//channel := make(chan integrationOutput)
		go evalIntegralN(Distance, contactR, var_list, times_list, numVariances, d_ind, channel)
		/*
		   // With this code structure, output is not set until the channel is ready to
		   // send. Hence we get no performance gain bc the code is still running
		   // entirely sequentially. Instead, we should only keep the go routine calls
		   // in this for loop, and then in a subsequent for loop, receive the results
		   // as they finish.
		   output := <- channel
		   index := output.d_index
		   CDF_arrays[index] = output.CDFs
		   Integrated_CDF_arrays[index] = output.integratedCDFs
		   //fmt.Println("Index is :",index)
		   //fmt.Println("channel output is :\n", output)
		*/

	}

	for i := 1; i <= numDistances; i++ {
		output := <-channel // this line will execute as soon as the channel has something to send
		index := output.d_index
		CDF_arrays[index] = output.CDFs
		Integrated_CDF_arrays[index] = output.integratedCDFs

		//fmt.Println("channel output is :\n", output)
		fmt.Println("d_ind ", index, " completed.")
	}
	close(channel)
	return CDF_arrays, Integrated_CDF_arrays
}

func evalIntegralN(Distance, contactR float64, var_list, times_list []float64, numVariances, d_ind int, Ch chan integrationOutput) {
	var Output integrationOutput
	CDF_vector := make([]float64, numVariances)
	Integrated_CDF_vector := make([]float64, numVariances)

	numTermsInApprox := 200 // Include *up to* this many terms
	numConsecutive := 20    //10
	//nan_reached := 0

	for v_ind := 0; v_ind < numVariances; v_ind++ {
		//for v_ind := 0; v_ind < 20; v_ind++{
		Variance := var_list[v_ind]
		theta := math.Pi * rand.Float64()
		phi := 2.0 * math.Pi * rand.Float64()
		muX := Distance * math.Sin(theta) * math.Cos(phi)
		muY := Distance * math.Sin(theta) * math.Sin(phi)
		muZ := Distance * math.Cos(theta)
		//fmt.Println("Variance/muX/muY/muZ = ",Variance,muX,muY,muZ)

		// Theorem 4.2b.1 Mathai-Provost p.95

		E := []float64{Variance, Variance, Variance} // diagonal of covariance matrix

		/*
		   % Steps:
		   0. compute b_vector
		   1. compute d(k) for k = 1:5
		   2. compute c(k-1)

		   Follow steps on p.96

		   Since the covariance matrix is diagonal,
		   P <- 3x3 Identity
		   D <- Diagonal matrix of eigenVals. All eigenVals equal to sigma^2_NEW
		   i.e. THE EIGENVALUES ARE JUST THE VARIANCES ALONG EACH DIMENSION AND IN OUR
		   ISOMETRIC CASE, ALL VARIANCES ARE EQUAL
		*/

		// b = mu / sqrt(variance) [ see 3.1a.3 in Mathai & Provost]
		b := []float64{(float64(1) / math.Sqrt(E[1])) * muX, (float64(1) / math.Sqrt(E[1])) * muY, (float64(1) / math.Sqrt(E[1])) * muZ}

		y := float64(contactR * contactR)

		chunks := make([]int, (numTermsInApprox / 10))
		chunks[0] = 1

		for j := 1; j < len(chunks); j++ {
			chunks[j] = chunks[j-1] + 10
		}
		//fmt.Println("chunks ",chunks)

		d := make([]float64, numTermsInApprox)
		c := make([]float64, numTermsInApprox)
		c[0] = math.Exp(-0.5*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])) * (1.0 / math.Sqrt(2*E[0])) * (1.0 / math.Sqrt(2*E[1])) * (1.0 / math.Sqrt(2*E[2]))
		//fmt.Println("c[0],b,E are: ",c[0],b,E)
		//fmt.Println("muX, muY, muZ, Distance are: ",muX, muY, muZ, Distance)
		Fs := make([]float64, numTermsInApprox)
		Fs[0] = c[0] * math.Pow(y, (float64(3.0)/float64(2.0))) / math.Gamma((float64(3.0)/float64(2.0))+1.0)

		limit_found := false
		converged := 0
		for mP := 0; mP < len(chunks)-1; mP++ {

			if limit_found == true {
				break
			}

			// Compute CDFs (i.e. Fs) for next chunk of k values and determine if
			// it has converged.
			for k := chunks[mP]; k < chunks[mP+1]; k++ {
				//fmt.Println("k = ",k)
				d[k] = 0.5 * ((1.0-float64(k)*(math.Pow(b[0], 2)))/math.Pow(2*E[0], float64(k)) + (1.0-float64(k)*(math.Pow(b[1], 2)))/math.Pow(2*E[1], float64(k)) + (1.0-float64(k)*(math.Pow(b[2], 2)))/math.Pow(2*E[2], float64(k)))

				// the formula for c(k) on p93 involves a sum from 0 to k-1
				for r := 0; r < k; r++ {
					c[k] = c[k] + d[k-r]*c[r]
				}
				c[k] = c[k] / float64(k)

				// The current approximation of the CDF to k terms
				// is a sum of k terms:
				Fs[k] = Fs[k-1] + math.Pow(-1.0, float64(k))*c[k]*math.Pow(y, (float64(3.0)/float64(2.0))+float64(k))/math.Gamma((float64(3.0)/float64(2.0))+float64(k)+1.0)
				//fmt.Println("F[k],k,c[k] are: ", Fs[k],k,c[k])

				if Fs[k] > 1 {
					panic("CDF value greater than 1 or less than zero")
				}

				// Test if CDF has converged. I.e., are last numConsecutive terms roughly equal?
				// I.e. are successive differences less than 1e-5?

				//fmt.Println("Fs[k] =",Fs[k])
				//fmt.Println("c[k] = ",c[k])
				//fmt.Println(math.Pow(-1.0, float64(k)) * c[k] * math.Pow(y, (float64(3.0)/float64(2.0)) + float64(k)) / math.Gamma((float64(3.0)/float64(2.0)) + float64(k) + 1.0))

				// UPDATE: DISREGRADE COMMENT PARAGRAPH BELOW.
				// Sometimes Fs[k] converges to a small magnitude negative value.
				// This happens with numConsecutive=20. Further, The c[k]
				// grow quite large in magnitude while the difference between
				// Fs[k-1] and Fs[k] (i.e. the long formula above) grows ever smaller.
				// This leads me to believe the convergence is genuine even though a
				// CDF should not be negative. Previously, I handled the possibly of a negative
				// limit by taking the max of (Fs[k],0). However this lead to all the
				// CDF[v_ind] elements being zero after some particular v_ind. The resulting
				// integrated CDF curve then plateaued. I am now making the choice
				// to take the absolute value of Fs[k] if it's actually converging
				// to avoid the possibility of negative CDF values.
				if math.Abs(Fs[k]-Fs[k-1]) < 1.0e-5 {
					converged++
					if converged >= numConsecutive {
						limit_found = true
						CDF_vector[v_ind] = math.Max(Fs[k], 0) //math.Abs(Fs[k])
						//fmt.Println("CDF_vector[v_ind]",CDF_vector[v_ind])
						//if v_ind >= 1 {
						//    fmt.Println("CDF ratio is: ", CDF_vector[v_ind-1]/CDF_vector[v_ind]) }
						//fmt.Println("Fs[k] =",Fs[k])

						// Trapezoidal Integration
						if v_ind > 0 {
							//fmt.Println(times_list[v_ind],times_list[v_ind-1])
							//fmt.Println("integral is ",(times_list[v_ind]-times_list[v_ind-1])*(CDF_arrays[d_ind][v_ind] + CDF_arrays[d_ind][v_ind-1])/2)
							Integrated_CDF_vector[v_ind] = Integrated_CDF_vector[v_ind-1] + (times_list[v_ind]-times_list[v_ind-1])*(CDF_vector[v_ind]+CDF_vector[v_ind-1])/2
							//fmt.Println("Integrated_CDF_vector[v_ind]",Integrated_CDF_vector[v_ind])

						}
					}
				} else {
					converged = 0
				}

				if limit_found {
					break
				}

			}

		}
	}

	// Gather the results into a struct
	// Send the struct over the channel, Ch
	// Don't close the chnnel here.
	Output.CDFs = CDF_vector
	Output.integratedCDFs = Integrated_CDF_vector
	Output.d_index = d_ind
	Ch <- Output

}

func writeCSV(data [][]float64, fileName string) bool {
	file, err := os.Create(fileName)
	checkError("Cannot create file", err)
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	for row, _ := range data {
		str := []string{}
		for _, element := range data[row] {
			text := fmt.Sprintf("%g", element)
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

func writeCSV2(data []float64, fileName string) bool {
	file, err := os.Create(fileName)
	checkError("Cannot create file", err)
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	str := []string{}
	for _, element := range data {
		text := fmt.Sprintf("%g", element)
		str = append(str, text)
	}

	err = writer.Write(str)
	checkError("Cannot write to file", err)

	return true
}

func MAXelem(f []float64) (int, float64) {
	max := 0.0
	index := 0
	for i, val := range f {
		if val > max {
			max = val
			index = i
		}
	}
	return index, max
}

func maxSecondDerivative(f []float64, h1, h2 float64, numTimePoints1, numTimePoints2 int) float64 {
	if numTimePoints1+numTimePoints2 != len(f) {
		panic("len(numTimePoints1) + len(numTimePoints2) should equal len(f)")
	}

	// The time point spacing used to evaluate f (i.e. the CDFs) is constant for
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
			DDf1[n] = (f1[n+1] - 2.0*f1[n] + f1[n-1]) / math.Pow(h1, 2)
		}
	}

	DDf2 := make([]float64, len(f2))
	for n := 1; n < len(f2)-1; n++ {
		if f2[n-1] == 0 {
			continue
		} else {
			DDf2[n] = (f2[n+1] - 2.0*f2[n] + f2[n-1]) / math.Pow(h2, 2)
		}
	}

	DDf := append(DDf1, DDf2...)
	_, maxDDf := MAXelem(DDf)
	return maxDDf
}
