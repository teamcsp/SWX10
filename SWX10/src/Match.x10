import x10.io.Console;
import x10.io.Reader;
import x10.io.File;
import x10.lang.Int;

class Match {
  public static def main(args:Rail[String]) {
	val SIZE:Long = 127;
	var HEAD_PARSED:Boolean = Boolean.FALSE;
	
    val acids:Rail[Char] = new Rail[Char](SIZE, '0'); 
    Console.OUT.printf("Acid char rail size is %d\n", acids.size);
    
    val inputMatrixFileName = "matrices/BLOSUM62";
    val inputMatrix = new File(inputMatrixFileName);
    val i = inputMatrix.openRead();
    for (line in i.lines()) {
      if ( line.trim().charAt(0 as Int) != '#' && line != null ) {
    	if (HEAD_PARSED == Boolean.FALSE){
    	  HEAD_PARSED = Boolean.TRUE;
    	}
    	parseLine(line.trim(), acids, HEAD_PARSED);
        Console.OUT.println("HEAD PARSED IS " + HEAD_PARSED);
      }
    }
    
    for (acid in acids) {
    	Console.OUT.println(acid + " ");
    }
  }

  private static def parseLine(line:String, acids:Rail[Char], HEAD_PARSED:Boolean):void {
    var splitedStrings:Rail[String] = line.split("  ");
    Console.OUT.println(splitedStrings);
    for (var i:Long=0;i < splitedStrings.size; i++)
    {
    	Console.OUT.print(i);
    	acids(i) = splitedStrings(i)(0 as Int);
    }
  }
}
