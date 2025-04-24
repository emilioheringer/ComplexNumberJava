﻿# JavaComplexLib

A Java library for handling complex numbers with full support for trigonometric and hyperbolic functions.

## ✨ Features

- Representation of complex numbers in **rectangular** and **polar** form.
- Basic operations: addition, subtraction, multiplication, division.
- Advanced math: exponential, logarithmic, powers.
- Trigonometric functions: `sin`, `cos`, `tan`, and their inverses.
- Hyperbolic functions: `sinh`, `cosh`, `tanh`, and more.

## 📦 Installation

You can download the JAR from the [Releases](https://github.com/emilioheringer/JavaComplexLib/releases) page or use a build system like **Gradle** or **Maven**.

### ✅ Gradle

```gradle
dependencyResolutionManagement {
    repositoriesMode.set(RepositoriesMode.FAIL_ON_PROJECT_REPOS)
    repositories {
        mavenCentral()
        maven { url 'https://jitpack.io' }
    }
}

dependencies {
    implementation 'com.github.emilioheringer:JavaComplexLib:1.0.0'
}
```

### ✅ Maven

```xml
<repositories>
    <repository>
        <id>jitpack.io</id>
        <url>https://jitpack.io</url>
    </repository>
</repositories>

<dependency>
    <groupId>com.github.emilioheringer</groupId>
    <artifactId>JavaComplexLib</artifactId>
    <version>main-d67721cb63-1</version>
</dependency>
```

## 🧪 Example Usage

```java
import com.heringer.Complex;

public class Main {
    public static void main(String[] args) {
        Complex z1 = new Complex(3, 4); // 3 + 4i
        Complex z2 = new Complex(1, -2); // 1 - 2i

        Complex sum = z1.add(z2);
        Complex product = z1.multiply(z2);

        sum.showRec()
        sum.showPolar()
    }
}
```

## 📚 Documentation

Full documentation is available at:  
👉 [https://emilioheringer.github.io/ComplexJavadoc/com/heringer/package-summary.html](https://emilioheringer.github.io/ComplexJavadoc/com/heringer/package-summary.html)

[![](https://jitpack.io/v/emilioheringer/ComplexNumberJava.svg)](https://jitpack.io/#emilioheringer/ComplexNumberJava)
