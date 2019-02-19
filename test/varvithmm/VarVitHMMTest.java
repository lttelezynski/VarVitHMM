/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package varvithmm;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author lucas
 */
public class VarVitHMMTest {
    VarVitHMM VarVitHMMInstance;
    
    public VarVitHMMTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
        System.out.println("VarVitHMMTest: Before method setUp()");
        VarVitHMMInstance = new VarVitHMM();
    }
    
    @After
    public void tearDown() {
        System.out.println("VarVitHMMTest: After method tearDown()");
        VarVitHMMInstance = null;
    }

    /**
     * Test of main method, of class VarVitHMM.
     */
//    @Test
//    public void testMain() {
//        System.out.println("main");
//        String[] args = null;
//        VarVitHMM.main(args);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

    /**
     * Test of VarVitHMM method, of class VarVitHMM.
     */
//    @Test
//    public void testVarVitHMM() {
//        System.out.println("VarVitHMM");
//        String[] args = null;
//        VarVitHMM instance = new VarVitHMM();
//        instance.VarVitHMM(args);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    
    @Test
    public void test_log_zero() {
        double expResult = -Double.MAX_VALUE;
        double result = VarVitHMMInstance.log(0);
        System.out.println("***VarVitHMMTest test_log_zero()");
        Assert.assertEquals(expResult, result, 1E300);
    }

    @Test
    public void test_log_one() {
        double expResult = 0;
        double result = VarVitHMMInstance.log(1);
        System.out.println("***VarVitHMMTest test_log_one()");
        
        Assert.assertEquals(expResult, result, 1E300);
    }

    @Test
    public void test_log_E() {
        double expResult = 1;
        double result = VarVitHMMInstance.log(Math.E);
        System.out.println("***VarVitHMMTest test_log_E()");
        
        Assert.assertEquals(expResult, result, 1E300);
    }
    

}
