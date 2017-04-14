#include "Variant.hpp"




// ============================================================================
void Variant::fromString (const std::string& rep)
{
    switch (type)
    {
        case 'i': intVal = std::stoi (rep); break;
        case 'd': doubleVal = std::stod (rep); break;
        case 's': stringVal = rep; break;
        case 'n': break;
    }
}

Variant::operator int() const
{
    switch (type)
    {
        case 'i': return intVal;
        case 'd': return doubleVal;
        case 's': return std::stoi (stringVal);
        default: return 0;
    }
}

Variant::operator double() const
{
    switch (type)
    {
        case 'i': return intVal;
        case 'd': return doubleVal;
        case 's': return std::stod (stringVal);
        default: return 0.0;
    }
}

Variant::operator std::string() const
{
    switch (type)
    {
        case 'i': return std::to_string (intVal);
        case 'd': return std::to_string (doubleVal);
        case 's': return stringVal;
        default: return "<null>";
    }
}

void Variant::printToStream (std::ostream& stream) const
{
    switch (type)
    {
        case 'i': stream << intVal; break;
        case 'd': stream << doubleVal; break;
        case 's': stream << stringVal; break;
        default: stream << "<null>"; break;
    }
}

std::ostream& operator<< (std::ostream& os, const Variant& var)
{
    var.printToStream (os);
    return os;
}

std::ostream& operator<< (std::ostream& os, const Variant::NamedValues& namedValues)
{
    for (auto& entry : namedValues)
    {
        os << entry.first << ": " << entry.second << std::endl;;
    }
    return os;
}

Variant::NamedValues Variant::fromCommandLine (int argc, const char* argv[])
{
    Variant::NamedValues namedValues;

    for (int n = 0; n < argc; ++n)
    {
        std::string arg = argv[n];
        std::string::size_type indexOfEqualSign = arg.find ('=');

        if (indexOfEqualSign != std::string::npos)
        {
            std::string key = arg.substr (0, indexOfEqualSign);
            std::string val = arg.substr (indexOfEqualSign + 1);
            namedValues[key] = val;
        }
    }
    return namedValues;
}

void Variant::update (Variant::NamedValues& target, const Variant::NamedValues& source)
{
    for (auto& entry : source)
    {
        if (target.find (entry.first) == target.end())
        {
            throw std::runtime_error ("unrecognized '" + entry.first + "'");
        }
        else
        {
            target[entry.first].fromString (entry.second);
        }
    }
}

void Variant::updateFromCommandLine (Variant::NamedValues& target, int argc, const char* argv[])
{
    Variant::NamedValues source = fromCommandLine (argc, argv);
    update (target, source);
}
