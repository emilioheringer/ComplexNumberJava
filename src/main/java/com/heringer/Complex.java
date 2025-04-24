package com.heringer;

/**
 * Represents a complex number with real and imaginary parts.
 * Supports both rectangular and polar forms.
 * Provides static and instance methods for arithmetic operations on complex
 * numbers.
 */
public class Complex implements Icomplex {

    private double Re, Img;
    String Mode = "Rec"; // "Rec" for rectangular, "Pol" for polar

    /**
     * Default constructor. Initializes a complex number as 0 + 0i.
     */
    public Complex() {
    }

    /**
     * Constructs a complex number using specified mode.
     *
     * @param Re   the real part or magnitude depending on the mode.
     * @param Img  the imaginary part or angle depending on the mode.
     * @param mode "Pol" for polar form, anything else for rectangular.
     */
    public Complex(Number Re, Number Img, String mode) {
        if (mode.toString().equals("Pol")) {
            double num_real = Re.doubleValue() * Math.cos(Img.doubleValue());
            double num_imag = Re.doubleValue() * Math.sin(Img.doubleValue());

            this.Img = num_imag;
            this.Re = num_real;
            this.Mode = "Pol";
        } else {
            this.Img = Img.doubleValue();
            this.Re = Re.doubleValue();
        }
    }

    /**
     * Constructs a complex number in rectangular form.
     *
     * @param Re  real part.
     * @param Img imaginary part.
     */
    public Complex(Number Re, Number Img) {
        this.Img = Img.doubleValue();
        this.Re = Re.doubleValue();
    }

    /**
     * Constructs a complex number from a string.
     * Supports both rectangular (e.g., "3+4i") and polar (e.g., "5<60") formats.
     *
     * @param complex string representation of a complex number.
     */
    public Complex(String complex) {
        complex = complex.replace("i", "").replace(" ", "");

        if (complex.contains("<")) {
            int angleIndex = complex.indexOf('<');
            String parteReal = complex.substring(0, angleIndex);
            String parteImaginaria = complex.substring(angleIndex + 1);

            double modulo = Double.parseDouble(parteReal);
            double anguloGraus = Double.parseDouble(parteImaginaria);
            double anguloRadianos = Math.toRadians(anguloGraus);

            this.Re = modulo * Math.cos(anguloRadianos);
            this.Img = modulo * Math.sin(anguloRadianos);
            this.Mode = "Pol";
        } else {
            int operadorIndex = Math.max(complex.indexOf('+'), complex.indexOf('-', 1));
            if (operadorIndex == -1) {
                throw new IllegalArgumentException("Invalid format. Use 'a + bi' or 'r<θ'");
            }

            char operador = complex.charAt(operadorIndex);
            String parteReal = complex.substring(0, operadorIndex);
            String parteImaginaria = complex.substring(operadorIndex + 1);

            this.Re = Double.parseDouble(parteReal);
            this.Img = Double.parseDouble(parteImaginaria);

            if (operador == '-') {
                this.Img = -this.Img;
            }

            this.Mode = "Rec";
        }
    }

    /**
     * Prints the complex number in rectangular form.
     */
    @Override
    public void showRec() {
        String sinal = Img >= 0 ? "+" : "-";
        System.out.printf("%.2f %s %.2fi%n", Re, sinal, Math.abs(Img));
    }

    /**
     * Prints the complex number in polar form (magnitude < angle).
     */
    @Override
    public void showPolar() {
        double z = Math.sqrt(Re * Re + Img * Img);
        double theta = Math.atan2(Img, Re);
        double thetaGraus = Math.toDegrees(theta);

        System.out.printf("%.2f < %.2f°%n", z, thetaGraus);
    }

    /** @return the real part of the complex number */
    public double getReal() {
        return this.Re;
    }

    /** @return the imaginary part of the complex number */
    public double getIma() {
        return this.Img;
    }

    /** @return the magnitude (modulus) of the complex number */
    public double getMagnetude() {
        return Math.sqrt(this.Re * this.Re + this.Img * this.Img);
    }

    /** @return the angle (in degrees) of the complex number */
    public double getAngleDegrees() {
        return Math.toDegrees(Math.atan2(this.Img, this.Re));
    }

    /** @return the angle (in radians) of the complex number */
    public double getAngleRad() {
        return Math.atan2(this.Img, this.Re);
    }

    /**
     * Sums multiple complex numbers given as strings.
     *
     * @param complex variable number of strings representing complex numbers.
     * @return resulting Complex after summation.
     */
    public static Complex sum(String... complex) {
        double re = 0;
        double img = 0;

        if (complex == null || complex.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        for (String complexStr : complex) {
            if (complexStr == null || complexStr.trim().isEmpty()) {
                throw new IllegalArgumentException("Invalid complex number string: " + complexStr);
            }

            Complex a = new Complex(complexStr);
            re += a.getReal();
            img += a.getIma();
        }
        return new Complex(re, img);
    }

    /**
     * Sums multiple complex numbers.
     *
     * @param complexNumbers variable number of Complex instances.
     * @return resulting Complex after summation.
     */
    public static Complex sum(Complex... complexNumbers) {
        double re = 0;
        double img = 0;

        if (complexNumbers == null || complexNumbers.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        for (Complex c : complexNumbers) {
            if (c == null) {
                throw new IllegalArgumentException("Null complex number encountered.");
            }

            re += c.getReal();
            img += c.getIma();
        }

        return new Complex(re, img);
    }

    /**
     * Subtracts multiple complex numbers given as strings.
     *
     * @param complex variable number of strings representing complex numbers.
     * @return resulting Complex after subtraction.
     */
    public static Complex sub(String... complex) {
        double re = 0;
        double img = 0;

        if (complex == null || complex.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        for (String complexStr : complex) {
            if (complexStr == null || complexStr.trim().isEmpty()) {
                throw new IllegalArgumentException("Invalid complex number string: " + complexStr);
            }

            Complex a = new Complex(complexStr);
            re -= a.getReal();
            img -= a.getIma();
        }
        return new Complex(re, img);
    }

    /**
     * Subtracts multiple complex numbers.
     *
     * @param complexNumbers variable number of Complex instances.
     * @return resulting Complex after subtraction.
     */
    public static Complex sub(Complex... complexNumbers) {
        double re = 0;
        double img = 0;

        if (complexNumbers == null || complexNumbers.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        for (Complex c : complexNumbers) {
            if (c == null) {
                throw new IllegalArgumentException("Null complex number encountered.");
            }

            re -= c.getReal();
            img -= c.getIma();
        }

        return new Complex(re, img);
    }

    /**
     * Multiplies multiple complex numbers given as strings.
     *
     * @param complex variable number of strings representing complex numbers.
     * @return resulting Complex after multiplication.
     */
    public static Complex Prod(String... complex) {
        double mag = 1;
        double ang = 0;

        if (complex == null || complex.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        for (String complexStr : complex) {
            if (complexStr == null || complexStr.trim().isEmpty()) {
                throw new IllegalArgumentException("Invalid complex number string: " + complexStr);
            }

            Complex a = new Complex(complexStr);
            mag *= a.getMagnetude();
            ang += a.getAngleDegrees();
        }
        return new Complex(mag, ang, "Pol");
    }

    /**
     * Multiplies multiple complex numbers.
     *
     * @param complexNumbers variable number of Complex instances.
     * @return resulting Complex after multiplication.
     */
    public static Complex Prod(Complex... complexNumbers) {
        if (complexNumbers == null || complexNumbers.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        double mag = 1;
        double ang = 0;

        for (Complex c : complexNumbers) {
            if (c == null) {
                throw new IllegalArgumentException("Null complex number encountered.");
            }

            mag *= c.getMagnetude();
            ang += c.getAngleDegrees();
        }

        return new Complex(mag, ang, "Pol");
    }

    /**
     * Divides a sequence of complex numbers represented as strings.
     * The division is performed in the order given: result = c1 / c2 / c3 / ...
     *
     * @param complex variable number of strings representing complex numbers.
     * @return resulting Complex after division.
     * @throws IllegalArgumentException if any complex number has a magnitude of
     *                                  zero.
     */
    public static Complex Div(String... complex) {
        if (complex == null || complex.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        Complex result = new Complex(complex[0]);

        for (int i = 1; i < complex.length; i++) {
            Complex a = new Complex(complex[i]);

            // Check for division by zero
            if (a.getMagnetude() == 0) {
                throw new IllegalArgumentException("Cannot divide by a complex number with zero magnitude.");
            }

            double mag = result.getMagnetude() / a.getMagnetude();
            double ang = result.getAngleDegrees() - a.getAngleDegrees();
            result = new Complex(mag, ang, "Pol");
        }

        return result;
    }

    /**
     * Divides a sequence of Complex numbers.
     * The division is performed in the order given: result = c1 / c2 / c3 / ...
     *
     * @param complexNumbers variable number of Complex instances.
     * @return resulting Complex after division.
     * @throws IllegalArgumentException if any complex number has a magnitude of
     *                                  zero.
     */
    public static Complex Div(Complex... complexNumbers) {
        if (complexNumbers == null || complexNumbers.length == 0) {
            throw new IllegalArgumentException("At least one complex number must be provided.");
        }

        Complex result = complexNumbers[0];

        // Check for division by zero in the first number
        if (result.getMagnetude() == 0) {
            throw new IllegalArgumentException("Cannot divide by a complex number with zero magnitude.");
        }

        for (int i = 1; i < complexNumbers.length; i++) {
            Complex c = complexNumbers[i];

            // Check for division by zero
            if (c.getMagnetude() == 0) {
                throw new IllegalArgumentException("Cannot divide by a complex number with zero magnitude.");
            }

            double mag = result.getMagnetude() / c.getMagnetude();
            double ang = result.getAngleDegrees() - c.getAngleDegrees();
            result = new Complex(mag, ang, "Pol");
        }

        return result;
    }

    /**
     * Adds another complex number to this one (in-place).
     *
     * @param a Complex number to add.
     * @return this Complex after addition.
     */
    @Override
    public Complex add(Complex a) {
        this.Re += a.getReal();
        this.Img += a.getIma();
        return this;
    }

    /**
     * Subtracts another complex number from this one (in-place).
     *
     * @param a Complex number to subtract.
     * @return this Complex after subtraction.
     */
    @Override
    public Complex sub(Complex a) {
        this.Re -= a.getReal();
        this.Img -= a.getIma();
        return this;
    }

    /**
     * Multiplies the current Complex number with another Complex number.
     * The multiplication is done using the polar representation.
     *
     * @param a The Complex number to multiply with the current Complex number.
     * @return The current Complex number instance after multiplication.
     */
    @Override
    public Complex Prod(Complex a) {
        // Multiplying magnitudes and adding angles in polar form
        double mag = this.getMagnetude() * a.getMagnetude();
        double ang = this.getAngleDegrees() + a.getAngleDegrees();

        // Updating the current instance with the result of the multiplication
        this.Re = mag * Math.cos(Math.toRadians(ang)); // Convert angle back to radians for the real part
        this.Img = mag * Math.sin(Math.toRadians(ang)); // Convert angle back to radians for the imaginary part

        return this;
    }

    /**
     * Divides the current Complex number by another Complex number.
     * Division is performed in polar form.
     *
     * @param a The Complex number to divide by.
     * @return The current Complex number instance after division.
     * @throws ArithmeticException If dividing by a zero complex number.
     */
    @Override
    public Complex Div(Complex a) {
        // Check for division by zero (both real and imaginary parts should not be zero)
        if (a.getReal() == 0 && a.getIma() == 0) {
            throw new ArithmeticException("Cannot divide by zero.");
        }

        // Dividing magnitudes and subtracting angles in polar form
        double mag = this.getMagnetude() / a.getMagnetude();
        double ang = this.getAngleDegrees() - a.getAngleDegrees();

        // Updating the current instance with the result of the division
        this.Re = mag * Math.cos(Math.toRadians(ang)); // Convert angle back to radians for the real part
        this.Img = mag * Math.sin(Math.toRadians(ang)); // Convert angle back to radians for the imaginary part

        return this;
    }

    /**
     * Conjugate Complex
     *
     * @param a complex number
     */
    public static Complex Conj(Complex a) {
        return new Complex(a.Re, -a.Img);
    }

    /**
     * Conjugate Complex of the instance of the object complex
     *
     */
    @Override
    public Complex getConj() {
        return new Complex(this.Re, -this.Img);
    }

    /**
     * Applies the Riemann projection to the current complex number.
     * This maps the complex number onto the Riemann sphere.
     *
     * @return A new Complex object representing the projected value.
     */
    @Override
    public Complex riemannProjection() {
        double modSquared = this.Re * this.Re + this.Img * this.Img + 1.0;
        double projectedRe = 2.0 * this.Re / modSquared;
        double projectedImg = 2.0 * this.Img / modSquared;
        return new Complex(projectedRe, projectedImg);
    }

    /**
     * Calculates the complex exponential of this number.
     * That is, returns e^(a + bi) = e^a * (cos(b) + i*sin(b)).
     *
     * @return A new Complex representing e^(this).
     */
    public Complex exp() {
        double expRe = Math.exp(this.Re);
        double realPart = expRe * Math.cos(this.Img);
        double imagPart = expRe * Math.sin(this.Img);
        return new Complex(realPart, imagPart);
    }

    /**
     * Returns the natural logarithm (ln) of this complex number.
     * The branch cut lies along the negative real axis (standard branch).
     *
     * @return A new Complex representing ln(z)
     */
    public Complex log() {
        double magnitude = Math.hypot(Re, Img); // mesmo que abs()
        double angle = Math.atan2(Img, Re); // mesmo que arg()
        return new Complex(Math.log(magnitude), angle);
    }

    /**
     * Computes the base-10 logarithm of this complex number.
     *
     * @return A new Complex representing log10(this)
     */
    public Complex log10() {
        Complex lnZ = this.log(); // Logaritmo natural
        double log10Real = lnZ.getReal() / Math.log(10);
        double log10Imag = lnZ.getIma() / Math.log(10);
        return new Complex(log10Real, log10Imag);
    }

    /**
     * Raises this complex number to the power of an integer exponent.
     *
     * @param n The integer exponent
     * @return A new Complex representing this^n
     */
    public Complex pow(int n) {
        double magnitude = Math.pow(this.getMagnetude(), n); // Magnitude elevada ao expoente
        double angle = this.getAngleDegrees() * n; // Multiplica o ângulo pelo expoente
        return new Complex(magnitude, angle, "Pol");
    }

    /**
     * Calculates the square root of this complex number.
     *
     * @return A new Complex representing the square root of this number
     */
    public Complex sqrt() {
        double magnitude = Math.sqrt(this.getMagnetude()); // Raiz quadrada da magnitude
        double angle = this.getAngleDegrees() / 2; // Divide o ângulo por 2
        return new Complex(magnitude, angle, "Pol");
    }

    /**
     * Calculates the sine of this complex number.
     *
     * @return A new Complex representing sin(this)
     */
    public Complex sin() {
        double realPart = Math.sin(this.Re) * Math.cosh(this.Img);
        double imagPart = Math.cos(this.Re) * Math.sinh(this.Img);
        return new Complex(realPart, imagPart);
    }

    /**
     * Calculates the cosine of this complex number.
     *
     * @return A new Complex representing cos(this)
     */
    public Complex cos() {
        double realPart = Math.cos(this.Re) * Math.cosh(this.Img);
        double imagPart = -Math.sin(this.Re) * Math.sinh(this.Img);
        return new Complex(realPart, imagPart);
    }

    /**
     * Calculates the tangent of this complex number.
     *
     * @return A new Complex representing tan(this)
     */
    public Complex tan() {
        Complex sinVal = this.sin();
        Complex cosVal = this.cos();
        return sinVal.Div(cosVal);
    }

    /**
     * Calculates the inverse sine (arcsin) of this complex number.
     *
     * @return A new Complex representing asin(this)
     */
    public Complex asin() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex z = this.Prod(i); // Multiplica por i (número imaginário)
        Complex result = z.add(new Complex(1, 0)); // Soma 1
        result = result.sqrt(); // Raiz quadrada
        Complex logTerm = result.add(this.add(i).sqrt()).log(); // Logaritmo da expressão
        return logTerm.Prod(i.negate()); // Multiplica pelo negativo de i
    }

    /**
     * Calculates the inverse cosine (arccos) of this complex number.
     *
     * @return A new Complex representing acos(this)
     */
    public Complex acos() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex result = this.add(i);
        result = result.sqrt();
        Complex logTerm = result.add(this.add(i).sqrt()).log(); // Logaritmo da expressão
        return logTerm.Prod(i.negate());
    }

    /**
     * Calculates the inverse tangent (arctan) of this complex number.
     *
     * @return A new Complex representing atan(this)
     */
    public Complex atan() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex z = this.Prod(i); // Multiplica por i
        Complex result = z.add(new Complex(1, 0)); // Soma 1
        result = result.sqrt(); // Raiz quadrada
        Complex logTerm = result.add(this.add(i).sqrt()).log(); // Logaritmo da expressão
        return logTerm.Prod(i.negate());
    }

    /**
     * Returns the negation of this complex number.
     *
     * @return A new Complex representing -this
     */
    public Complex negate() {
        return new Complex(-this.Re, -this.Img);
    }

    /**
     * Calculates the hyperbolic sine of this complex number.
     *
     * @return A new Complex representing sinh(this)
     */
    public Complex sinh() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex expPlus = this.add(i).exp(); // e^(z+i)
        Complex expMinus = this.sub(i).exp(); // e^(z-i)

        return expPlus.sub(expMinus).Div(new Complex(2, 0)); // (e^(z+i) - e^(z-i)) / 2
    }

    /**
     * Calculates the hyperbolic cosine of this complex number.
     *
     * @return A new Complex representing cosh(this)
     */
    public Complex cosh() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex expPlus = this.add(i).exp(); // e^(z+i)
        Complex expMinus = this.sub(i).exp(); // e^(z-i)

        return expPlus.add(expMinus).Div(new Complex(2, 0)); // (e^(z+i) + e^(z-i)) / 2
    }

    /**
     * Calculates the hyperbolic tangent of this complex number.
     *
     * @return A new Complex representing tanh(this)
     */
    public Complex tanh() {
        Complex sinhValue = this.sinh();
        Complex coshValue = this.cosh();

        return sinhValue.Div(coshValue); // sinh(z) / cosh(z)
    }

    /**
     * Calculates the inverse hyperbolic sine of this complex number.
     *
     * @return A new Complex representing asinh(this)
     */
    public Complex asinh() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex term = this.add(this.Prod(i)); // z + i*z
        term = term.sqrt(); // sqrt(z + i*z)

        // ln(z + sqrt(z + i*z))
        return term.add(this).log();
    }

    /**
     * Calculates the inverse hyperbolic cosine of this complex number.
     *
     * @return A new Complex representing acosh(this)
     */
    public Complex acosh() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex term = this.add(this.Prod(i)); // z + i*z
        term = term.sqrt(); // sqrt(z + i*z)

        // ln(z + sqrt(z + i*z))
        return term.add(this).log();
    }

    /**
     * Calculates the inverse hyperbolic tangent of this complex number.
     *
     * @return A new Complex representing atanh(this)
     */
    public Complex atanh() {
        Complex i = new Complex(0, 1); // Unidade imaginária
        Complex term = this.add(i); // z + i
        Complex term2 = this.sub(i); // z - i

        // 0.5 * ln((z + i) / (z - i))
        return term.Div(term2).log().Prod(new Complex(0.5, 0));
    }

}

