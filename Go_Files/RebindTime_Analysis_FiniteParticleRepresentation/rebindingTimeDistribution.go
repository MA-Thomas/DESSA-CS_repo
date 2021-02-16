package main

import (
	"fmt"
	"math"
	"math/rand"
)

type rebindOutput struct {
	waitT float64
}

func rebindDistrib(IntegratedPropensity, times []float64) ([]float64, []float64, int) {
	// Each bin will count the number of wait times in the corresponding time range
	// I.e. Bins[t] will count wait times between BinTimes[t] and BinTimes[t+1]
	numBins := int(2e8)//int(1e7)
	Bins := make([]float64, numBins)
	numTimes := len(times)

	// Define the bin times used to assign wait times to bins
	// Bins are logarithmically spaced. I.e. each successive entry is larger than
	// the prev entry by the same factor.
	// UPDATE: TO COMPUTE PDF FROM HISTOGRAM, WE NEED BINS OF EQUAL WIDTH. BINWIDTH IS USED IN THE CALCULATION.
	//         THEREFORE, DO NOT USE LOGARITHMICALLY SPACED BINS.
	//factor := math.Pow((times[len(times)-1]/times[0]), 1.0/float64(numBins))
	constant := (times[len(times)-1] - times[0]) / float64(numBins)
	BinTimes := make([]float64, numBins)
	BinTimes[0] = times[0]
	for t := 0; t < numBins-1; t++ {
		//BinTimes[t+1] = BinTimes[t]*factor
		BinTimes[t+1] = BinTimes[t] + constant
	}




	// Go runs out of memory when numTrials is too big.
	// Use a lower numTrials and just repeat multiple times.

	numRepeats := 1e7
	numTrials := 1e3
	//waitTimes := make([][]float64,int(numRepeats) )
	var waitTime float64


	for z := 0; z < int(numRepeats); z++ {



		//waitTimes[z] = make([]float64, int(numTrials))



		// Have Go workers compute the wait times
		channel := make(chan rebindOutput, int(numTrials))
		for i := 0; i < int(numTrials); i++ {

			go decide(numTimes, IntegratedPropensity, times, channel)
		}

		// Collect results from Go workers
		for i := 0; i < int(numTrials); i++ {
			output := <-channel // this line will execute as soon as the channel has something to send
			//waitTimes[z][i] = output.waitT
			waitTime = output.waitT

			// Use binarySearch() to quickly find the right bin for the wait time.
			//ans, index := binarySearch(waitTimes[z][i], BinTimes)
			ans, index := binarySearch(waitTime, BinTimes)
			if ans {
				// waitTimes[i] is assigned to the appropriate bin.
				Bins[index]++
			} else {
				// waitTimes[i] does not fit in any bin. It is larger than the last time
				// we are considering, e.g. 1 or 2 seconds.
				// In this case, I'm deciding to increment the last bin.
				Bins[len(Bins)-1]++
				//fmt.Println("not found")
			}
			/*
			if waitTimes[i] <= BinTimes[len(BinTimes)-1] {


				for n := 0; n < int(numBins)-1; n++ {
					if waitTimes[i] > BinTimes[n] && waitTimes[i] <= BinTimes[n+1] {
						Bins[n]++
						break
					}
				}

			}
			*/
		}
	}
	//return Bins, BinTimes, waitTimes, numBins
	return Bins, BinTimes, numBins
}

func decide(numTimes int, IntegratedPropensity, times []float64, Ch chan rebindOutput) {
	var Output rebindOutput
	Pk := math.Log(1.0 / rand.Float64())

	w := -1.0 // Default value

	// Use binarySearch() to quickly find the wait time.
	ans, index := binarySearch(Pk, IntegratedPropensity)
	if ans {
		w = times[index]
	}
/*
	for t := 0; t < numTimes; t++ {
		//fmt.Println(i, t)
		if IntegratedPropensity[t] >= Pk {

			w = times[t]
			//fmt.Println(waitTimes[i])
			break
		}
	}
*/

	// Gather the results into a struct
	// Send the struct over the channel, Ch
	// Don't close the chnnel here.
	Output.waitT = w
	Ch <- Output

}


// Binary Search.
// Returns the (bool,int). If bool==true, then int is the index of the first
// element of haystack to exceed threshold.
func binarySearch(threshold float64, haystack []float64) (bool, int) {


	low := 0
	high := len(haystack) - 1

	if haystack[low] == threshold {
		return true, low
	}
	if haystack[high] == threshold {
		return true, high
	}

	if haystack[high] < threshold {
		return false, 666
	}

	if haystack[low] > threshold {
		return false, 666
	}

	for low <= high{
		median := (low + high) / 2

		if haystack[median] < threshold {

			if haystack[median+1] >= threshold {
				return true, median+1
			} else {

			low = median + 1
			}
		}else{
			high = median - 1
		}
	}

	//if low == len(haystack) || haystack[low] != threshold{
	//	return false, 666
	//}

	return true, low
}
