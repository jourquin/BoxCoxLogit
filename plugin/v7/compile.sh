#!/bin/sh
# Compile the sample modal-split plugins and generate the jar files

# Add the Nodus main jar and libs to the classpath (EDIT ACCORDING TO YOUR PATH)
export NODUS7_HOME=/Applications/Nodus7

# Compile the source code of the plugins
javac -cp $NODUS7_HOME/nodus7.jar:$NODUS7_HOME/lib/* -source 1.8 -target 1.8 BoxCox1.java
javac -cp $NODUS7_HOME/nodus7.jar:$NODUS7_HOME/lib/* -source 1.8 -target 1.8 BoxCox2.java 

# Create the JAR files
jar cf BoxCox1.jar BoxCox1.class
jar cf BoxCox2.jar BoxCox2.class

# Remove the compiled class file (not mandatory)
rm BoxCox1.class
rm BoxCox2.class
