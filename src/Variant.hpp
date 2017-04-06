#ifndef Variant_hpp
#define Variant_hpp

#include <string>
#include <map>
#include <iostream>


class Variant
{
public:
    typedef std::map<std::string, Variant> NamedValues;

    Variant() : type ('n') {}
    Variant (int val) : intVal (val), type ('i') {}
    Variant (double val) : doubleVal (val), type ('d') {}
    Variant (std::string val) : stringVal (val), type ('s') {}
    Variant (const char* val) : stringVal (val), type ('s') {}

    void fromString (const std::string& rep);
    void printToStream (std::ostream& stream) const;
    operator int() const;
    operator double() const;
    operator std::string() const;

    static NamedValues fromCommandLine (int argc, const char* argv[]);
    static void updateFromCommandLine (Variant::NamedValues& target, int argc, const char* argv[]);
    static void update (Variant::NamedValues& target, const Variant::NamedValues& source);

private:
    int intVal;
    double doubleVal;
    std::string stringVal;
    char type;
};

std::ostream& operator<< (std::ostream& os, const Variant& var);
std::ostream& operator<< (std::ostream& os, const Variant::NamedValues& namedValues);

#endif // Variant_hpp
