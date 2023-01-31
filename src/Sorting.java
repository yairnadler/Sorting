import java.util.Random;

import Plotter.Plotter;


public class Sorting{

	final static int SELECT_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int COUNTING_VS_QUICK_LENGTH = 16;
	final static int BUBBLE_VS_MERGE_LENGTH = 12;
	final static int MERGE_VS_QUICK_SORTED_LENGTH =11;
	final static double T = 600.0;
	
	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 * 
	 * @param arr - the array to be sorted
	 */
    public static void quickSort(double[] arr){
        quickSort(arr, 0, arr.length - 1);
    }

    public static void quickSort(double[] arr, int p, int r){
        if (p < r - 2){
            int q = partition(arr, p, r);
            quickSort(arr, p, q-1);
            quickSort(arr, q+1, r);
        }
        else{
            bubbleSort(arr, p ,r);
        }
        
    }
    
    public static int partition(double[] arr, int p, int r){
        double pivot = arr[r];
        int i = p - 1;
		for (int j = p; j <= r - 1; j++){
			if (arr[j] <= pivot){
				i++;
				swap(arr, i, j);
			}
		}
		swap(arr, i+1, r);
		return (i+1);
    }
	


	/**
	 * Given an array arr and an index i returns the the i'th order statstics in arr. 
	 * In other words, it returns the element with rank i in the array arr.
	 * 
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * 
	 * Should run in average complexity of O(n), and worst case complexity of O(n^2)
	 * 
	 **/
    public static double QuickSelect(double[] arr, int i){
		double result = QuickSelect(arr, 0, arr.length - 1, i);
        return result;
	}

    public static double QuickSelect(double[] arr, int p, int r, int i){
    	if (i == 0) i++;
        if (p == r) return arr[p];
        int q = partition(arr, p, r);
        int m = q - p + 1;
        if (i == m) return arr[q];
        if (i < m) return QuickSelect(arr, p, q - 1, i);
        else return QuickSelect(arr, q+1, r, i - m);
    }
	
	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
    public static void mergeSort(double[] arr){
        if (arr.length == 0){
            throw new Error("Empty array!");
        }
        megreSort(arr, 0, arr.length - 1);    
    }

    public static void megreSort(double[] arr, int p, int r){
        if (p < r){
            int q = (p + r) / 2;
            megreSort(arr, p, q);
            megreSort(arr, q + 1, r);
            merge(arr, p, q, r);
        }
    }
    /**
     * merges two sorted subarrays into array
     * @param arr - array
     * @param p - start index of left subarray
     * @param q - end index of left array and start index of right array
     * @param r - end index of right subarray
     */
    public static void merge(double[] arr, int p, int q, int r){
        int length1 = q - p + 1;
        int length2 = r - q;
        double[] leftArr = new double[length1+1];
        double[] rightArr = new double[length2+1];
        for (int i = 0; i < length1; i++){
            leftArr[i] = arr[i + p];
        }
        for (int i = 0; i < length2; i++){
            rightArr[i] = arr[i + q + 1];
        }
        leftArr[length1] = Double.POSITIVE_INFINITY;
        rightArr[length2] = Double.POSITIVE_INFINITY;
        int i = 0;
        int j = 0;
        for (int k = p; k <= r; k++){
            if (leftArr[i] <= rightArr[j]){
                arr[k] = leftArr[i];
                i++;
            }
            else{
                arr[k] = rightArr[j];
                j++;
            }
        }

    }

	
	/**
	 * Sorts a given array using bubble sort.
	 * 
	 * The algorithm should run in complexity O(n^2).
	 * 
	 * @param arr - the array to be sorted
	 */
    public static void bubbleSort(double[] arr){
        if (arr.length == 0){
            throw new Error("Empty array!");
        }
        int length = arr.length;
        int stop = 0;
        for (int i = 0; i < length; i++){
            int swap = 0;
            if (stop == 1) break;
            for (int j = 1; j < length - i; j++){
                if ((arr[j]) < arr[j-1]){
                    swap(arr, j-1, j);
                    swap = 1;
                }
            }
            if (swap == 0){
                stop = 1;
            }
        }
    }
    /**
     * same as bubble sort but for a portion of the given array
     * @param arr - given array
     * @param p - start index
     * @param r - end index
     */
    public static void bubbleSort(double[] arr, int p, int r){
        if (arr.length == 1) return;
        int stop = 0;
        for (int i = p; i <= r + 1; i++){
            int swap = 0;
            if (stop == 1) break;
            for (int j = p + 1; j < r + 1; j++){
                if ((arr[j]) < arr[j-1]){
                    swap(arr, j-1, j);
                    swap = 1;
                }
            }
            if (swap == 0){
                stop = 1;
            }
        }
    }
    

	/**
	 * Sorts a given array, using the counting sort algorithm.
	 * You may assume that all elements in the array are between 0 and k (not including k).
	 * 
	 * Should run in complexity O(n + k) in the worst case.
	 * 
	 * @param arr - an array with positive integers
	 * @param k - an upper bound for the values of all elements in the array.
	 */
    public static void countingSort(int[] arr, int k){
		int[] output = new int[arr.length];
        int[] temp = new int[k+1];
        for (int i = 1; i < arr.length; i++){
            temp[arr[i]]++;
        }
        for (int i = 1; i < temp.length; i++){
            temp[i] += temp[i-1];
        }
        for (int i = arr.length - 1; i >= 1; i--){
            output[temp[arr[i]]] = arr[i];
            temp[arr[i]]--;
        }
        for (int i = 0; i < arr.length; i++) {
        	arr[i] = output[i];
        }
        
	}
    
    /**
 	 * swaps between 2 elments in a given int array
 	 * @param arr - given array
 	 * @param i - first index
 	 * @param j - second index
 	 */
     public static void swap(double[] arr, int i, int j){
 		double temp = arr[i];
 		arr[i] = arr[j];
 		arr[j] = temp;
 	}
     
 	/**
 	 * swaps between 2 elments in a given int array
 	 * @param arr - given array
 	 * @param i - first index
 	 * @param j - second index
 	 */
 	public static void swap(int[] arr, int i, int j){
 		int temp = arr[i];
 		arr[i] = arr[j];
 		arr[j] = temp;
 	}

	/**
	 * prints a given array
	 * used for debugging
	 * @param arr = given array
	 */
 	
	 public static void printArr(double[] arr){
		 StringBuilder s = new StringBuilder();
		 s.append("[");
		 for (int i = 0; i < arr.length; i++){
			 s.append(arr[i] + ", ");
		 }
		 s.deleteCharAt(s.length() - 1);
		 s.deleteCharAt(s.length() - 1);
		 s.append("]");
		 System.out.println(s.toString());
	 }

	/**
	 * prints a given array
	 * used for debugging
	 * @param arr = given array
	 */
	 public static void printArr(int[] arr){
		StringBuilder s = new StringBuilder();
		s.append("[");
		for (int i = 0; i < arr.length; i++){
			s.append(arr[i] + ", ");
		}
		s.deleteCharAt(s.length() - 1);
		s.deleteCharAt(s.length() - 1);
		s.append("]");
		System.out.println(s.toString());
	}
	/**
	 * This method runs the sorting algorithms on random arrays and
	 * presents a graph of the running times.
	 * To run a method, uncomment the call to it in the main method.
	 */
	public static void main(String[] args) {
	
		//countingVsQuick();
		//mergeVsQuick();
		//mergeVsQuickOnSortedArray();
		//mergeVsBubble();
		//QuickSelectVsQuickSort();
		
	}
	


	private static void countingVsQuick() {
		double[] quickTimes = new double[COUNTING_VS_QUICK_LENGTH];
		double[] countingTimes = new double[COUNTING_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < COUNTING_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumCounting = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				int[] b = new int[size];
				for (int j = 0; j < a.length; j++) {
					b[j] = r.nextInt(size);
					a[j] = b[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				countingSort(b, size);
				endTime = System.currentTimeMillis();
				sumCounting += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			countingTimes[i] = sumCounting/T;
		}
		Plotter.plot("Counting sort on arrays with elements < n", countingTimes, "Quick sort on arrays with elements < n", quickTimes);
		
	}


	
	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick(){
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on random array", quickTimes, "merge sort on random array", mergeTimes);
	}
	
	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray(){
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge  += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}
	/**
	 * Compares merge sort and bubble sort on random arrays
	 */
	public static void mergeVsBubble(){
		double[] mergeTimes = new double[BUBBLE_VS_MERGE_LENGTH];
		double[] bubbleTimes = new double[BUBBLE_VS_MERGE_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < BUBBLE_VS_MERGE_LENGTH; i++) {
			long sumMerge = 0;
			long sumBubble = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				mergeSort(a);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
				startTime = System.currentTimeMillis();
				bubbleSort(b);
				endTime = System.currentTimeMillis();
				sumBubble += endTime - startTime;
			}
			mergeTimes[i] = sumMerge/T;
			bubbleTimes[i] = sumBubble/T;
		}
		Plotter.plot("merge sort on random array", mergeTimes, "bubble sort on random array", bubbleTimes);
	}
	
	



	/**
	 * Compares the quick select algorithm with a random rank, and the quick sort algorithm.
	 */
	public static void QuickSelectVsQuickSort(){
		double[] QsortTimes = new double[SELECT_VS_QUICK_LENGTH];
		double[] QselectTimes = new double[SELECT_VS_QUICK_LENGTH];
		Random r = new Random();
		long startTime, endTime;
		for (int i = 0; i < SELECT_VS_QUICK_LENGTH; i++) {
			long sumQsort = 0;
			long sumQselect = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQsort += endTime - startTime;
				startTime = System.currentTimeMillis();
				//printArr(b);
				QuickSelect(b, r.nextInt(size));
				endTime = System.currentTimeMillis();
				sumQselect += endTime - startTime;
			}
			QsortTimes[i] = sumQsort/T;
			QselectTimes[i] = sumQselect/T;
		}
		Plotter.plot("quick sort with an arbitrary pivot", QsortTimes, "quick select with an arbitrary pivot, and a random rank", QselectTimes);
	}
	

}
