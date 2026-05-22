#ifndef JABBER_INTERPOLANT
#define JABBER_INTERPOLANT

#include <map>
#include <vector>
#include <span>

namespace jabber
{

/**
 * @brief Base interface for a simple R->R continuous function.
 */
class Function
{
protected:
public:

   /// Evaluate function at \p x.
   virtual double operator() (const double &x) const = 0;
   virtual ~Function() = default;
};

/// Piecewise-linear interpolant.
class PWLinear : public Function
{
private:

   /// Struct of y-value and slope associated with previous PW line.
   struct Line
   {
      double y;
      double m;
   };

   /// Map of Key = x, Value=Line (y-value for x, and slope).
   std::map<double, Line> func_map_;

protected:

   /// Get const reference to \ref func_map_.
   const std::map<double, Line>& Map() const { return func_map_; }

public:
   /// Construct piecewise linear interpolant through ( \p x_k , \p y_k ).
   PWLinear(std::span<const double> x_k, std::span<const double> y_k);

   double operator() (const double &x) const override;
};

/// Piecewise log-log interpolant (linear on log-log scaling).
class PWLogLog : public Function
{
private:

   /// Struct of y-value and log-space slope associated with previous PW line.
   struct Line
   {
      double y;
      double m;
   };

   /// Map of Key = x, Value=Line (y-value for x, and log-space slope).
   std::map<double, Line> func_map_;

protected:

   /// Get const reference to \ref func_map_.
   const std::map<double, Line>& Map() const { return func_map_; }

public:
   /**
    * @brief Construct piecewise log-log (linear in log10 space) interpolant
    * through ( \p x_k , \p y_k ).
    */
   PWLogLog(std::span<const double> x_k, std::span<const double> y_k);
   
   double operator() (const double &x) const override;
};

} // namespace jabber

#endif // JABBER_INTERPOLANT
