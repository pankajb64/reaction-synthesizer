/**
 * @author pankaj
 *
 */

/**
 * This class ReactionSynthesizer is a model class which approximates
 * possible chemical reactions generated during various computing constructs.
 * The following constructs are currently supported - 
 * 
 * NOT - See abs() and abs2step()
 * Increment - See increment()
 * Decrement - See decrement()
 * Copy - See copy()
 * If a > n - See ifGreaterThan()
 * Selection Sort - See sort()
 * 
 * Based on the paper "Rate-Independent Constructs for Chemical Computation" by "Phillip Senum, Marc Riedel, 
 * Department of Electrical and Computer Engineering, University of Minnesota, Minneapolis, Minnesota, United States of America"
 * Citation: Senum P, Riedel M (2011) Rate-Independent Constructs for Chemical Computation. 
 * PLoS ONE 6(6): e21414. doi:10.1371/journal.pone.0021414
 * 
 * @author pankaj
 *
 */
public class ReactionSynthesizer
{

	private int g;
	private int t;
	private int f;
	
	public ReactionSynthesizer()
	{
		g = 0;
		t = 0;
		f = 0;
	}
	
	/**
	 * Implements Selection Sort Algorithm. Please read comments for more description.
	 * @param length
	 */
	public void sort(int length)
	{
		String i = "i";
		String n = "n";
		String n1 = "n(c)";
		String n2 = "n(c)(c)";
		String j = "j";
		String min = "min";
		String Amin = "A(min)";
		String Aj = "A(j)";
		String AminAb = abs( Amin, true );
		String AjAb = abs( Aj, true );
		String phi = phi();
		String minAb = abs( min, true);
		String jc = "jc(c)";
		String t = "";
		String gP = "gP", gPP = "gPP", gPP1 = "gPP1", gPP2 = "gPP2";
		String g1, g2, g11, g12, g13, g14, g15;
		String jPAb = abs( prime( j ), true );
		String minPAb = abs( prime(min), true);
		String temp = "temp";
		String Ai = "A(i)";
		String jAb = abs( j, true );
		String iPAb = abs( prime( i ), true );
		String F = "F";
		
		copy(n1, n);
		decrement( n1 );
		copy(i, n1);
		
		
		//outer while loop starts - while i > 0
		copy(n2, i);
		decrement(n2);
		copy(j, n2);
		copy(min, i);
		
		
		//inner while loop starts - while j > 0
		
		t = IfGreaterThan( Amin, Aj );
		
		//These three reactions are for assigning min to j if Amin > Aj. They destroy min and then copy jc to min
		copy(jc, j); 
		printBiMolecular( t, AjAb, min, AjAb, t, false );
		printBiMolecular( t, minAb, AjAb, minAb, t, false );
		copy(min, jc);
		
		decrement( j );
		
		//Generators for inner while loop
		g1 = "g"+(g-1);
		g2 = "g"+(g-2);
		
		printTriMolecular( j, jPAb, minPAb, j, gP, false );
		printBiMolecular( gP, jc, gP, false );
		printUniMolecular( gP, g1, g2, false ); 
		
		//End of Inner while loop
		
		
		/*
		 * We can't really synthesize reactions for swapping elements, it has no meaning.
		 * The value of min at the end of this iteration should be used to as an index
		 * to pick the element at that position. We have added an extra fuel "F" below
		 * which will be used to clear the value of min at the end of outer loop iteration.
		 * The outer loop will remain paused in a sense until the fuel is added everytime.
		 */
		decrement( i );
		
		//Generators for outer while loop
		g11="g"+(g-7);
		g12="g"+(g-6);
		g13="g"+(g-5);
		g14="g"+(g-4);
		g15="g"+(g-3);
		
		/* 
		 * we generate signal gPP here and clear min and n2 using the fuel F.
		 * Note that there is no need to clear j as it will be already 0 at the end
		 * of inner while loop.
		 */
		printTriMolecular( i, jAb, iPAb, i, gPP, false );
		printTriMolecular( F, gPP, min, F, gPP, false );
		printTriMolecular( minAb, n2 , gPP, minAb, gPP, false );
		printUniMolecular( gPP, g11, gPP1, false );
		printBiMolecular( gPP1, g11, g12, gPP2, false );
		printBiMolecular( gPP2, g12, g13, g14, g15, false );

		System.out.println(  );
	}
	
	/**
	 * Implements the comparator - If a is greater than b
	 * @param a The first parameter
	 * @param b The second parameter
	 * @return The truth signal indicator t
	 */
	public String IfGreaterThan(String a, String b)
	{
		String ac = a+"(c)";
		String bc = b+"(c)";
		String phi = phi();
		String fuel = nextFuel();
		String t = nextTruth();
		
		copy(ac, a);
		copy(bc, b);
		
		String acAb = abs(ac, true);
		String bcAb = abs( bc, true );
		
		printBiMolecular( ac, bc, phi, true );
		printTriMolecular( fuel, ac, bcAb, ac, bcAb, t, false );
		printTriMolecular( acAb, bc, t, acAb, bc, false );
		printTriMolecular( acAb, bcAb, t, acAb, bcAb, false );
		
		System.out.println( );
		
		return t;
	}
	
	/**
	 * Implements the copier - Copies "from" into "to"
	 * @param to The element where from should be copied to 
	 * @param from The element to be copied
	 */
	public void copy(String to, String from)
	{
		String g = nextSignal();
		String fromP = prime(from);
		String fromAb = abs(from, true);
		String gAb = abs(g, true);
		String phi = phi();
		
		printBiMolecular( g, from, g, fromP, false );
		printBiMolecular( g, fromAb, phi, false );
		printBiMolecular( gAb, fromP, from, to, false );
		
		System.out.println(  );
	}
	
	/**
	 * Common Unity for Increment and Decrement
	 * @param s The string to be incremented or decremented
	 * @param increment true if increment, false otherwise
	 * @return The signal used to increment or decrement
	 */
	public String incDec(String s, boolean increment)
	{
		String g = nextSignal();
		String sP = prime( s );
		String sP2 = prefix( 2, sP );
		String sab = abs(s, true);
		String gpAb = abs2Step( g, true );
		String phi = phi();
		String sr = rx(s);
		String srab = abs(sr, true);
		String s2 = prefix( 2, s );
		
		printBiMolecular( s, g, sP, g, false );
		printBiMolecular( g, sab, phi, false );
		printBiMolecular( gpAb, sP2, s, sP, sr, true );
		printUniMolecular( sr, phi, false );
		
		if(increment)
		{
			printTriMolecular(srab, sP, gpAb, s2, false);
		}
		else
		{
			printTriMolecular(srab, sP, gpAb, phi, false);
		}
		
		System.out.println(  );
		
		return g;
	}
	
	/**
	 * Implements the incrementer - Increments the element by 1
	 * @param s The element to be incremented
	 * @return The signal used to increment
	 */
	public String increment(String s)
	{
		return incDec( s, true );
		
	}
	
	/**
	 * Implements the decrementer - Decrements the element by 1
	 * @param s The element to be decremented
	 * @return The signal used to decrement
	 */
	public String  decrement(String s)
	{
		return incDec( s, false );
	}
	
	/**
	 * Utility method to get the next signal notation
	 * @return the next signal
	 */
	public String nextSignal()
	{
		return "g"+g++;
	}
	
	/**
	 * Utility method to get the next fuel notation
	 * @return the next fuel
	 */
	public String nextFuel()
	{
		return "f"+f++;
	}
	
	/**
	 * Utility method to return the next truth notation
	 * @return the next truth
	 */
	public String nextTruth()
	{
		return "t"+t++;
	}
	
	/**
	 * Utility method to return the prime notation for the given element
	 * @param s the element whose prime should be generated
	 * @return the generated prime
	 */
	public String prime(String s)
	{
		return s+"'";
	}
	
	/**
	 * Implements the NOT gate - Creates an absence indicator for the given element
	 * @param s The string whose obscene indicator is to be generated
	 * @param printReactions true if the reactions should be printed, false otherwise
	 * @return The generated absence indicator
	 */
	public String abs(String s, boolean printReactions)
	{
		String abs = s+"(ab)"; 
		

		if(printReactions)
		{
			String phi = phi();
			String abs2 = prefix(2, abs);

			printUniMolecular( phi, abs, false );
			printBiMolecular( s, abs, s, true );
			printUniMolecular( abs2, abs, true );
			System.out.println(  );
		}
		
		return abs;
	}
	
	/**
	 * Implements the NOT gate - a 2 step method to generate the absence indicator for the given element.
	 * @param s the element whose absence indicator is to be generated
	 * @param printReactions true if the reactions should be printed, false otherwise
	 * @return the generated absence indicator
	 */
	public String abs2Step(String s, boolean printReactions)
	{
		String abs = abs( s, printReactions );
		String absP = abs(prime( s ), false);
		
		if(printReactions)
		{
			String abs2P = prefix(2, absP);

			printUniMolecular( abs, absP, false );
			printBiMolecular( s, absP, s, true );
			printUniMolecular( abs2P, absP, true );
			System.out.println(  );
		}
		
		return absP;
	}
	
	/**
	 * Utility method to generate the phi notation
	 * @return the phi notaton
	 */
	public String phi()
	{
		return "0";
	}
	
	/**
	 * Utility to prefix an element with a number
	 * @param prefix the number to be prefixed
	 * @param s the element to attach the prefix to
	 * @return the prefixed element
	 */
	public String prefix(int prefix, String s)
	{
		return prefix+s;
	}
	
	/**
	 * Utilty method to generate the rx of element
	 * @param s the element
	 * @return the rx notation
	 */
	public String rx(String s)
	{
		return s+"(rx)";
	}
	
	/**
	 * Utility method to print reactions with a single reactant and two products
	 * @param reactant the reactant
	 * @param product1 the first product
	 * @param product2 the second product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printUniMolecular(String reactant, String product1, String product2, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf( "\t\t\t\t\t\t\t\t\t%-20s\t------>\t%20s\t+\t%20s\t\t\t\t\t(%s)\n", reactant, product1, product2, rate );		
	}
	
	/**
	 * Utility method to print reactions with a single reactant and a single products
	 * @param reactant the reactant
	 * @param product the product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printUniMolecular(String reactant, String product, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf( "\t\t\t\t\t\t\t\t\t%-20s\t------>\t%20s\t\t\t\t\t\t\t\t\t(%s)\n", reactant, product, rate );
	}
	
	/**
	 * Utility method to print reactions with two reactants and three products
	 * @param reactant1 the first reactant
	 * @param reactant2 the second reactant
	 * @param product1 the first product
	 * @param product2 the second product
	 * @param product3 the third product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printBiMolecular(String reactant1, String reactant2, String product1, String product2, String product3, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf("\t\t\t\t\t%-20s\t+\t%-20s\t------>\t%20s\t+\t%20s\t+\t%20s\t(%s)\n", reactant1, reactant2, product1, product2, product3, rate);
	}
	
	/**
	 * Utility method to print reactions with two reactants and two products
	 * @param reactant1 the first reactant
	 * @param reactant2 the second reactant
	 * @param product1 the first product
	 * @param product2 the second product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printBiMolecular(String reactant1, String reactant2, String product1, String product2, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf("\t\t\t\t\t%-20s\t+\t%-20s\t------>\t%20s\t+\t%20s\t\t\t\t\t(%s)\n", reactant1, reactant2, product1, product2, rate);
	}
	
	/**
	 * Utility method to print reactions with two reactants and a single product
	 * @param reactant1 the first reactant
	 * @param reactant2 the second reactant
	 * @param product the product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printBiMolecular(String reactant1, String reactant2, String product, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf("\t\t\t\t\t%-20s\t+\t%-20s\t------>\t%20s\t\t\t\t\t\t\t\t\t(%s)\n", reactant1, reactant2, product, rate );
	}
	
	/**
	 * Utility method to print reactions with three reactants and a single product
	 * @param reactant1 the first reactant
	 * @param reactant2 the second reactant
	 * @param reactant3 the third reactant
	 * @param product the product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printTriMolecular(String reactant1, String reactant2, String reactant3, String product, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf("\t%-20s\t+\t%-20s\t+\t%20s\t------>\t%20s\t\t\t\t\t\t\t\t\t(%s)\n", reactant1, reactant2, reactant3, product, rate);
	}
	
	/**
	 * Utility method to print reactions with three reactants and two products
	 * @param reactant1 the first reactant
	 * @param reactant2 the second reactant
	 * @param reactant3 the third reactant
	 * @param product1 the first product
	 * @param product2 the second product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printTriMolecular(String reactant1, String reactant2, String reactant3, String product1, String product2, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf("\t%-20s\t+\t%-20s\t+\t%20s\t------>\t%20s\t+\t%20s\t\t\t\t\t(%s)\n", reactant1, reactant2, reactant3, product1, product2, rate);
	}
	
	/**
	 * Utility method to print reactions with three reactants and three products
	 * @param reactant1 the first reactant
	 * @param reactant2 the second reactant
	 * @param reactant3 the third reactant
	 * @param product1 the first product
	 * @param product2 the second product
	 * @param product3 the third product
	 * @param fast rate indicating whether reaction is fast or slow
	 */
	public void printTriMolecular(String reactant1, String reactant2, String reactant3, String product1, String product2, String product3, boolean fast)
	{
		String rate =  fast == true ? "fast" : "slow";
		System.out.printf("\t%-20s\t+\t%-20s\t+\t%20s\t------>\t%20s\t+\t%20s\t+\t%20s\t(%s)\n", reactant1, reactant2, reactant3, product1, product2, product3, rate);
	}
	
	/**
	 * The main method
	 * @param args System arguments
	 */
	public static void main(String[] args)
	{
		ReactionSynthesizer rs = new ReactionSynthesizer();
		//s.copy( "a", "b" );
		//s.increment( "x" );
		//s.decrement( "y" );
		//s.IfGreaterThan( "a", "b" );
		System.out.println( "Reactions to sort an array of size 10\n\n\n" );
		rs.sort( 10 );
	}
}
