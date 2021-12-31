package demo;

/**
 * Created by pdoviet on 9/9/2015.
 */
class A implements Cloneable {
    public int i = 10;
}

class B extends A {
    public int i = 20;

    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }
}

public class Test {
    public static void main(String[] args) throws CloneNotSupportedException {
        B b1 = new B();
        B b2 = (B) b1.clone();
        System.out.println(b1 == b2);
    }
}