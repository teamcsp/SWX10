import x10.io.Console;
import x10.io.Reader;
import x10.io.File;
import x10.lang.Int;
import x10.array.Array_2;
import x10.util.Timer;
import x10.util.RailBuilder;
import x10.lang.Math;
import x10.util.Pair;

class Match {
	public static def main(args:Rail[String]) {
		val SIZE:Long = 127;
		
		/*
		 * Arguments:
		 * --------------------------------------------
		 * fastaOne: fasta file containing the first sequence
		 * fastaTwo: fasta file containing the second sequence
		 * (fasta files are found in ./fasta directory)
		 * 
		 * blosum: 	file containing the Substitution Matrix
		 * (blosum files are found in ./matrices directory)
		 * 
		 * gapA:	gap opening penalty
		 * gapB:	gap extention penalty
		 */
		
		val fastaOne:String, fastaTwo:String, blosum:String, gapA:String, gapB:String;
		
		if(args.size == 5){
			fastaOne = args(0);
			fastaTwo = args(1);
			blosum = args(2);
			gapA = args(3);
			gapB = args(4);
		}
		else {
			// Some default parameters
			fastaOne = "AF043946.1";
			fastaTwo = "AF043947.1";
			blosum = "BLOSUM62";
			gapA = "10";
			gapB = "5";
		}
		
		// 1D Rail to store the raw (in chars) Substitution Matrix
		val acids:Rail[Char] = new Rail[Char](SIZE, '0'); 
		Console.OUT.printf("Acid char rail size is %d\n", acids.size);
		
		// 2D Array for Substitution Matrix
		val subMatrix: Array_2[Int] = new Array_2[Int](SIZE, SIZE);
		
		val inputMatrixFileName = "matrices/" + blosum;
		
		// Parsing the Matrix file
		val startTime: Long = Timer.nanoTime();
		parseMatrixFile(inputMatrixFileName, acids, subMatrix );
		val endTime: Long = Timer.nanoTime();
		
		Console.OUT.println("Time elapsed loading matrix : " + (endTime-startTime)/1000 + " microsecs");
		
		Console.OUT.print("Acid rail non 0 values : ");
		for (acid in acids) {
			if (acid != '0')
				Console.OUT.print(acid + " ");
		}
		Console.OUT.println();
		Console.OUT.println("Test sub matrix for subMatrix[M][B] expects -3; From matrix = " + subMatrix('M'.ord(),'B'.ord()));
		Console.OUT.println("Test sub matrix for subMatrix[B][B] expects 4; From matrix = " + subMatrix('B'.ord(),'B'.ord()));
		Console.OUT.println("Test sub matrix for subMatrix[Z][G] expects -2; From matrix = " + subMatrix('Z'.ord(),'G'.ord()));
		
		// Parse fasta files into stringA and stringB
		val stringA: Rail[Char];
		val stringB: Rail[Char];
		
		stringA = parseFasta(fastaOne);
		// for (char in stringA){
		// 	Console.OUT.print(char + " ");
		// }
		// Console.OUT.println();
		stringB = parseFasta(fastaTwo);
		// for (char in stringB){
		// 	Console.OUT.print(char + " ");
		// }
		Console.OUT.println();    
		
		// smith-waterman sequential version init params
		val scoringMatrix: Array_2[Long] = new Array_2[Long](stringA.size+1, stringB.size+1);
		val gapPenalty:Long = Long.parse(gapA) + Long.parse(gapB);
		var globalMax:Long = 0;
		var max_i:Long = -1, max_j:Long = -1;
		val parent2D: Array_2[Pair[Long,Long]] = new Array_2[Pair[Long, Long]](stringA.size+1, stringB.size+1);
		
		Console.OUT.println("Scoring Matrix Rows: " + scoringMatrix.numElems_1);
		Console.OUT.println("Scoring Matrix Cols: " + scoringMatrix.numElems_2);
		
		val result: Pair[Long, Pair[Long,Long]] = SmithWatermanSeq(subMatrix,stringA, stringB, scoringMatrix, gapPenalty, globalMax, max_i, max_j, parent2D);
		
		globalMax = result.first;
		max_i = result.second.first;
		max_j = result.second.second;
		
		Console.OUT.println("Max value in scoring matrix: " + globalMax);
		Console.OUT.println("Max i: " + max_i + " Max j: " + max_j);
		
		// Traceback from the element with the biggest score in the Scoring Matrix
		Console.OUT.println("Parent2D [max_i][max_j]: " + parent2D(max_i,max_j).first + " " + parent2D(max_i,max_j).second);
		
		var matchA:String = ""; var matchB:String = "";
		
		matchA += stringA(max_i-1);
		matchB += stringB(max_j-1);
		
		var tempI:Long = parent2D(max_i, max_j).first; 
		var tempJ:Long = parent2D(max_i, max_j).second;
		
		// Gap Count for counting the number of gaps encountered during traceback
		var gapCount: Long = 0;
		
		// Trace Length for tracking the length of matchA and matchB during traceback
		// Trace Length starts at 1, since the first char is already added to the matchA and matchB strings
		var traceLength:Long = 1; 
		
		// Identity Count for counting the true char matches between matchA and matchB
		// Identity Count starts at 1, assuming the first char from the two strings are a true match
		// i.e stringA(max_i-1) == stringB(max_j-1)
		var identityCount:Long = 1;
		
		while (parent2D(tempI, tempJ) != Pair(0,0)) {
			
			// Console.OUT.println("Current parent2D : " + parent2D(tempI, tempJ));
			
			val parent: Pair[Long, Long] = parent2D(tempI, tempJ);
			
			// Console.OUT.println("parent score: " + scoringMatrix(tempI,tempJ));
			
			// if parent is diagonal no problem
			if ( (tempI-1)== parent.first as Long && (tempJ-1) == parent.second as Long){
				matchA = stringA(tempI-1) + matchA;
				matchB = stringB(tempJ-1) + matchB;
				
				// Important to check if the chars are a true match
				if(stringA(tempI-1) == stringB(tempJ-1)){
					identityCount++;
				}
			}
			else {
				// parent is left
				if (parent.first as Long == tempI && parent.second as Long == (tempJ-1)) {
					// MatchB not affected
					matchB = stringB(tempJ-1) + matchB;
					
					// Add '-' for gap in MatchA
					matchA = "-" + matchA;
				} 
				// parent is top
				else {
					//MatchA not affected
					matchA = stringA(tempI-1) + matchA;
					
					// Add '-' for gap in MatchB
					matchB = "-" + matchB;
				}
				
				gapCount++;
			}
			
			traceLength++;
			
			tempI = parent.first as Long;
			tempJ = parent.second as Long;
			// Console.OUT.println("tempI : " + tempI + " tempJ : " + tempJ);
		}
		
		Console.OUT.println("Identity : " + identityCount + "/" + traceLength);
		Console.OUT.println("Gaps : " + gapCount + "/" + traceLength);
		Console.OUT.println("Match A : " + matchA);
		Console.OUT.println("Match B : " + matchB);
		
		// TODO: SMITH-WATERMAN IN PARALLEL
		
		val scoringMatrixP : Array_2[Long] = new Array_2[Long](stringA.size+1, stringB.size+1, (i:long, j:long)=>i==0? 0: (j==0) ? 0: -1);
		
		Console.OUT.println("scoreMat(0,0): " + scoringMatrixP(0,0));
		Console.OUT.println("scoreMat(1,0): " + scoringMatrixP(1,0));
		Console.OUT.println("scoreMat(3,3): " + scoringMatrixP(3,3));
		Console.OUT.println("scoreMat(0,1): " + scoringMatrixP(0,1));
		
		val gapPenaltyP:Long = Long.parse(gapA) + Long.parse(gapB);
		var globalMaxP:Long = 0;
		var max_iP:Long = -1, max_jP:Long = -1;
		val parent2DP: Array_2[Pair[Long,Long]] = new Array_2[Pair[Long, Long]](stringA.size+1, stringB.size+1);
		// Done Matrix to track the progress of the computation
		val doneMatrix: Array_2[Long] = new Array_2[Long](stringA.size+1, stringB.size+1);
		
		// Change top row of doneMatrix to 1s
		//for(j:Long in 0..(doneMatrix.numElems_2 - 1)){
		//	doneMatrix(1, j) = 1;
		//}
		
		// Change first column of doneMatrix to 1s
		//for(i:Long in 0..(doneMatrix.numElems_1 - 1)){
		//	doneMatrix(i, 1) = 1;
		//}
		
		var numThreads:Long = stringA.size;
		
		// Launch the threads. Number of threads required is equal to the length of stringA
		finish { for(i in 1..(numThreads)) async {
			// i implies the row that the thread is operating on
			
			var stop:Boolean = false;
			var j:Long = 1;
			var numCols:Long = scoringMatrix.numElems_2 -1 ;
			var top_done: long = 0;
			var diag_done: long = 0;
			var side_done: long = 0;
			
			// Each thread iterates leftwards 
			while(j <= numCols){
				// Wait until the 3 required elements in the scoring matrix is ready
				//while(true){
				
				when(scoringMatrixP(i-1,j) != -1 && scoringMatrixP(i-1,j-1) != -1 && scoringMatrixP(i,j-1)!= -1){
					Console.OUT.println("i = " + i + "  j = " + j );
					val matchP:Long = scoringMatrixP(i-1, j-1) + subMatrix(stringA(i-1).ord(), stringB(j-1).ord());
					val sideGapP:Long = scoringMatrixP(i, j-1) - gapPenaltyP;
					val topGapP:Long = scoringMatrixP(i-1, j) - gapPenaltyP;
					var maxP:Long = 0;
					
					// Compute scoringMatrix(i,j)
					if (matchP > sideGapP) {
						if (matchP > topGapP) {
							if ( matchP > 0 ) {
								// match is the greatest and is positive
								maxP = matchP;
								//parent2D(i,j) = Pair(i-1, j-1);
							}
							else {
								maxP = 0;
							}
						} else {
							if (topGapP > 0) {
								// topGap is the greatest and is positive
								maxP = topGapP;
								//parent2D(i,j) = Pair(i-1, j);
							}
							else {
								maxP = 0;
							}
						}
					} else if (sideGapP > topGapP) {
						if (sideGapP > 0) {
							// sideGap is the greatest and is positive
							maxP = sideGapP;
							//parent2D(i,j) = Pair(i, j-1);
						}
						else { 
							maxP = 0;
						}
					} else {
						if (topGapP > 0) {
							// topGap is the greatest and is positive
							maxP = topGapP;
							//parent2D(i,j) = Pair(i-1, j);
						}
						else { 
							maxP = 0;
						}
					} 
					atomic {
						if (maxP >= globalMaxP) {
							globalMaxP = maxP;
						}
					}
					// Escape from this while loop
					//break;
				}
				//}
				
				j++;
			}
			
		}
		}
		Console.OUT.println("Global Max P is :" + globalMaxP);
		
		return;
		
	}
	
	public static def SmithWatermanSeq(subMatrix: Array_2[Int],stringA: Rail[Char], stringB: Rail[Char],scoringMatrix: Array_2[Long], gapPenalty: Long, globalMax:Long, max_i:Long, max_j:Long, parent2D:Array_2[Pair[Long,Long]]): Pair[Long, Pair[Long,Long]] {
		var result:Pair[Long, Pair[Long,Long]] = Pair(0 as Long, Pair(0 as Long, 0 as Long));
		for(var i:Long = 1; i<scoringMatrix.numElems_1; i++){
			for(var j:Long = 1; j<scoringMatrix.numElems_2; j++){
				
				// Compute the 3 required values
				val match:Long = scoringMatrix(i-1, j-1) + subMatrix(stringA(i-1).ord(), stringB(j-1).ord());
				val sideGap:Long = scoringMatrix(i, j-1) - gapPenalty;
				val topGap:Long = scoringMatrix(i-1, j) - gapPenalty;
				
				// Find the Max value
				var max:Long = 0;
				
				if (match > sideGap) {
					if (match > topGap) {
						if ( match > 0 ) {
							// match is the greatest and is positive
							max = match;
							parent2D(i,j) = Pair(i-1, j-1);
						}
						else {
							max = 0;
						}
					} else {
						if (topGap > 0) {
							// topGap is the greatest and is positive
							max = topGap;
							parent2D(i,j) = Pair(i-1, j);
						}
						else {
							max = 0;
						}
					}
				} else if (sideGap > topGap) {
					if (sideGap > 0) {
						// sideGap is the greatest and is positive
						max = sideGap;
						parent2D(i,j) = Pair(i, j-1);
					}
					else { 
						max = 0;
					}
				} else {
					if (topGap > 0) {
						// topGap is the greatest and is positive
						max = topGap;
						parent2D(i,j) = Pair(i-1, j);
					}
					else { 
						max = 0;
					}
				}
				
				// Assign the Max value to the scoring matrix
				scoringMatrix(i, j) = max;
				
				// Note down the global max value of the scoring matrix
				if(max >= globalMax){
					
					result = Pair(max, Pair(i,j));
					
				}
				
				// Console.OUT.print(max + "\t");
			}
			// Console.OUT.println("");
		}
		return result;
	}
	
	// Parses a given FASTA into a Rail of Chars
	public static def parseFasta(fastaName:String): Rail[Char] {
		val fastaFile = new File("fasta/" + fastaName + ".fasta");
		val fastaLines = fastaFile.lines();
		fastaLines.next(); // Skip the first line
		
		// Rail Builder for building the rail of chars
		val railBuilder: RailBuilder[Char] = new RailBuilder[Char]();
		
		// Add each char from each line to the railBuilder
		for(line in fastaLines){
			val trimmedLine = line.trim();
			for(var i:Int = (0 as Int); i<trimmedLine.length(); i++ ){
				railBuilder.add(line.charAt(i));	
			}
		}
		
		// Return the rail of chars
		return railBuilder.result();
	}
	
	// Parses file with inputMatrixFileName, and stores the result into acids and subMatrix
	// Designed for reading the substituition matrix file
	private static def parseMatrixFile(inputMatrixFileName: String, acids:Rail[Char], subMatrix: Array_2[Int]): void {
		val inputMatrix = new File(inputMatrixFileName);
		
		val i = inputMatrix.openRead();
		var HEAD_PARSED:Boolean = Boolean.FALSE;
		
		for (line in i.lines()) {
			
			// Skip first few lines (i.e. lines starting with '#')
			if ( line.trim().charAt(0 as Int) != '#' && line != null ) {
				parseLine(line.trim(), acids, HEAD_PARSED, subMatrix);
				if (HEAD_PARSED == Boolean.FALSE){
					HEAD_PARSED = Boolean.TRUE;
				}
				//Console.OUT.println("HEAD PARSED IS " + HEAD_PARSED);
			}
		}	    
	}
	
	// Parses line. Results is stored in acids and subMatrix
	private static def parseLine(line:String, acids:Rail[Char], HEAD_PARSED:Boolean, subMatrix: Array_2[Int]):void {
		var splitedStrings:Rail[String]; 
		
		// Spliting the string using spaces
		if (HEAD_PARSED == Boolean.FALSE) {
			// Header elements are seperated by double-spaces
			// i.e "A  R  N  D  .. "
			splitedStrings = line.split("  ");
		} else {
			splitedStrings = line.split(" ");
		}
		
		
		var acid: Char = 0 as Char; 
		var extraCount: Int  = 0 as Int;
		
		// Processing each element of the splited string
		for (var i:Long=0;i < splitedStrings.size; i++)
		{
			if (HEAD_PARSED == Boolean.FALSE) {
				acids(i) = splitedStrings(i)(0 as Int);
			} else {
				if (splitedStrings(i).length() == 0 as Int) {
					extraCount++;
					continue;
				}
				if (i == 0) {
					acid = splitedStrings(i)(0 as Int);
					// Console.OUT.println("Current Head is : " + acid + " ");
				}
				else {
					val j: Int = i as Int - extraCount;
					//Console.OUT.print( splitedStrings(i) + "[ acids(j-1) is " + acids(j-1)+ "]");
					if (acids(j-1) != '0') {
						//Console.OUT.println("submatrix[ " + acid + "(" + acid.ord() + ")][" + acids(j-1) + "(" + acids(j-1).ord() +")] is : " + Int.parse(splitedStrings(i)) + " ");
						subMatrix(acid.ord(), acids(j-1).ord()) = Int.parse(splitedStrings(i));
					}
				}
			}
		}
	}
}
