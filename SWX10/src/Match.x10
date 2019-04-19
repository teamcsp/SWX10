import x10.io.Console;
import x10.io.Reader;
import x10.io.File;
import x10.lang.Int;
import x10.array.Array_2;

class Match {
  public static def main(args:Rail[String]) {
	val SIZE:Long = 127;
	var HEAD_PARSED:Boolean = Boolean.FALSE;
	
    val acids:Rail[Char] = new Rail[Char](SIZE, '0'); 
    Console.OUT.printf("Acid char rail size is %d\n", acids.size);
    
    val subMatrix: Array_2[Int] = new Array_2[Int](SIZE, SIZE);
    
    val inputMatrixFileName = "matrices/BLOSUM62";
    val inputMatrix = new File(inputMatrixFileName);
    val i = inputMatrix.openRead();
    for (line in i.lines()) {
      if ( line.trim().charAt(0 as Int) != '#' && line != null ) {
    	parseLine(line.trim(), acids, HEAD_PARSED, subMatrix);
    	if (HEAD_PARSED == Boolean.FALSE){
    		HEAD_PARSED = Boolean.TRUE;
    	}
        Console.OUT.println("HEAD PARSED IS " + HEAD_PARSED);
      }
    }
    
    for (acid in acids) {
    	if (acid != '0')
    	Console.OUT.println(acid + " ");
    }
  }

  private static def parseLine(line:String, acids:Rail[Char], HEAD_PARSED:Boolean, subMatrix: Array_2[Int]):void {
	var splitedStrings:Rail[String]; 
	if (HEAD_PARSED == Boolean.FALSE) {
		splitedStrings = line.split("  ");
	} else {
		splitedStrings = line.split(" ");
	}
	
    Console.OUT.println(splitedStrings);
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
    			Console.OUT.println("Current Head is : " + acid + " ");
    		}
    		else {
    			val j: Int = i as Int - extraCount;
    			Console.OUT.print( splitedStrings(i) + "[ acids(j-1) is " + acids(j-1)+ "]");
    			if (acids(j-1) != '0') {
    				Console.OUT.println("submatrix[ " + acid + "(" + acid.ord() + ")][" + acids(j-1) + "(" + acids(j-1).ord() +")] is : " + Int.parse(splitedStrings(i)) + " ");
    				subMatrix(acid.ord(), acids(j-1).ord()) = Int.parse(splitedStrings(i));
    			}
    		}
    	}
    }
  }
}
