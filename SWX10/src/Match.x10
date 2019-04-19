import x10.io.Console;
import x10.io.Reader;
import x10.io.File;
import x10.lang.Int;
import x10.array.Array_2;
import x10.util.Timer;

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
    
    val inputMatrixFileName = "matrices/BLOSUM62";
    
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
	if (HEAD_PARSED == Boolean.FALSE) {
		splitedStrings = line.split("  ");
	} else {
		splitedStrings = line.split(" ");
	}
	
    //Console.OUT.println(splitedStrings);
    var acid: Char = 0 as Char; 
    var extraCount: Int  = 0 as Int;
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
    			//Console.OUT.println("Current Head is : " + acid + " ");
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
