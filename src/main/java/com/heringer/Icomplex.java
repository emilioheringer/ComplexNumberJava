package com.heringer;

public interface Icomplex {

    public void showRec();
    public void showPolar();
    public double getReal();
    public double getIma();
    public double getMagnetude();
    public double getAngleDegrees();
    public double getAngleRad();
    public Complex add(Complex a);
    public Complex sub(Complex a);
    public Complex Prod(Complex a);
    public Complex Div(Complex a);
    public Complex getConj();
    public Complex riemannProjection();

}

