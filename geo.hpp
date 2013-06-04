/******************************************************************************
 *                           Spatial Index Utilities
 *                  Data structures - Geometric primitives.
 *****************************************************************************/

#ifndef _GEO__H_
# define _GEO__H_

# include <limits>
# include <vector>
# include <deque>
# include <algorithm>
# include <cstdlib>
# include <cmath>
# include <iostream>
# include <functional>

namespace Geo {

  /// Lenient comparison function.
  template <typename T>
  inline
  bool isGeometricallyEqual(T a, T b) {
    return a == b;
  }
  template <>
  inline
  bool isGeometricallyEqual<double>(double a, double b) {
    return std::abs(a - b) < 1e-15;
  }

  /*************************************************
   *
   * Point data type, representation, coercion.
   *
   *************************************************/
  typedef long CompactPoint;

  /// Cross compact that minimize space.
  inline
  CompactPoint crossCompactPoint(int x, int y) {
    const unsigned long xL = x;
    const unsigned long yL = y;
    return
      (xL & 0x0000000fUL)     |
      ((yL & 0x0000000fUL) << 4) |
      ((xL & 0x000000f0UL) << 4) |
      ((yL & 0x000000f0UL) << 8) |
      ((xL & 0x00000f00UL) << 8) |
      ((yL & 0x00000f00UL) << 12) |
      ((xL & 0x0000f000UL) << 12) |
      ((yL & 0x0000f000UL) << 16) |
      ((xL & 0x000f0000UL) << 16) |
      ((yL & 0x000f0000UL) << 20) |
      ((xL & 0x00f00000UL) << 20) |
      ((yL & 0x00f00000UL) << 24) |
      ((xL & 0x0f000000UL) << 24) |
      ((yL & 0x0f000000UL) << 28) |
      ((xL & 0xf0000000UL) << 28) |
      ((yL & 0xf0000000UL) << 32);
  }

  inline
  void crossDecompactPoint(CompactPoint point, int& x, int& y) {
    x = (point & 0x000000000000000fLL) |
      ((point & 0x0000000000000f00LL) >> 4) |
      ((point & 0x00000000000f0000LL) >> 8) |
      ((point & 0x000000000f000000LL) >> 12) |
      ((point & 0x0000000f00000000LL) >> 16) |
      ((point & 0x00000f0000000000LL) >> 20) |
      ((point & 0x000f000000000000LL) >> 24) |
      ((point & 0x0f00000000000000LL) >> 28);

    y = ((point & 0x00000000000000f0LL) >> 4) |
      ((point & 0x000000000000f000LL) >> 8) |
      ((point & 0x0000000000f00000LL) >> 12) |
      ((point & 0x00000000f0000000LL) >> 16) |
      ((point & 0x000000f000000000LL) >> 20) |
      ((point & 0x0000f00000000000LL) >> 24) |
      ((point & 0x00f0000000000000LL) >> 28) |
      ((point & 0xf000000000000000LL) >> 32);
  }

  /// Points.
  template <typename CoordT>
  struct Point {
    typedef CoordT coord_type;

    Point() {}

    Point(coord_type x, coord_type y) {
      coords[0] = x; coords[1] = y;
    }

    Point(coord_type c[2]) {
      coords[0] = c[0]; coords[1] = c[1];
    }

    Point(CompactPoint p) {
      crossDecompactPoint(p, coords[0], coords[1]);
    }

    template <typename>
    friend std::ostream& operator<<(std::ostream& out, const Point& p);

    bool operator==(const Point &other) const {
      return isGeometricallyEqual<coord_type>(x(), other.x()) && isGeometricallyEqual<coord_type>(y(), other.y());
    }

    bool operator!=(const Point &other) const {
      return !isGeometricallyEqual<coord_type>(x(), other.x()) || !isGeometricallyEqual<coord_type>(y(), other.y());
    }

    coord_type x() const { return coords[0]; }
    coord_type y() const { return coords[1]; }

    coord_type coords[2];
  };

  template <typename coord_type>
  std::ostream& operator<<(std::ostream& out, const Point<coord_type>& p) {
    out << p.x() << "," << p.y();
    return out;
  }


  /// Counter-clockwise turn.
  ///
  /// Three points are a counter-clockwise turn if > 0, clockwise if
  /// < 0, and collinear if ccw = 0 because ccw is a determinant that
  /// gives the signed area of the triangle formed by a, b and c.
  template <typename T>
  inline
  double counterClockwiseTurn(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
    // don't integer-overflow
    double ax = a.x();
    double ay = a.y();
    return (b.x() - ax) * (c.y() - ay) - (b.y() - ay) * (c.x() - ax);
  }

  /// Collinearity function.
  template <typename T>
  inline
  bool isCollinear(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
    return isGeometricallyEqual(counterClockwiseTurn(a, b, c), 0.0);
  }

  /// Dot product
  /// Can integer-overflow if the points are far away from each other!
  template <typename T>
  inline
  T dot_product(const Point<T>& a, const Point<T>&b, const Point<T>& c) {
    return (b.x() - a.x()) * (c.x() - b.x()) + (b.y() - a.y()) * (c.y() - b.y());
  }

  /*************************************************
   *
   * Vectors.
   *
   *************************************************/
  template <typename Coord>
  class Segment {
  public:
    Segment(const Segment& other)
    : _from(other._from), _to(other._to) {}
    Segment(const Point<Coord>& from, const Point<Coord>& to)
    : _from(from), _to(to) {}
    Segment(): _from(0, 0), _to(0, 0) {};

    const Point<Coord>& from() const {
      return _from;
    }
    const Point<Coord>& to() const {
      return _to;
    }

    bool operator==(const Segment& other) const {
      return _from == other._from
        && _to == other._to;
    }

  private:
    Point<Coord> _from;
    Point<Coord> _to;
  };


  /*************************************************
   * Regions.
   *
   * Disks.
   * Polygons, as:
   * - a closed set of Vectors.
   * - a set of Points.
   * Bounding Boxes
   * Union of Regions.
   * Intersection of Regions.
   *
   *************************************************/

  /// Axis-aligned Bounding Box.
  template <typename coord_type>
  class BoundingBox : public std::pair<Point<coord_type>,Point<coord_type> > {
  public:
    typedef Point<coord_type>                point_type;
    typedef std::pair<point_type,point_type> base_type;
    typedef BoundingBox<coord_type>          bounding_box_type;

    BoundingBox() {}

    BoundingBox(coord_type xMin, coord_type yMin, coord_type xMax, coord_type yMax)
    : base_type(point_type(xMin, yMin), point_type(xMax, yMax)) {}

    BoundingBox(coord_type min[2], coord_type max[2])
    : base_type(point_type(min[0], min[1]), point_type(max[0], max[1])) {}

    BoundingBox(const point_type& low, const point_type& high)
    : base_type(std::make_pair(low, high)) {}

    coord_type xMin() const { return this->first.x(); }
    coord_type yMin() const { return this->first.y(); }
    coord_type xMax() const { return this->second.x(); }
    coord_type yMax() const { return this->second.y(); }

  private:
    template <typename> friend class Polygon;
  };


  template <typename coord_type>
  class Disk {
  public:
    typedef Point<coord_type>       point_type;
    typedef BoundingBox<coord_type> bounding_box_type;

    Disk(point_type c, double r)
    : center(c), radius(r) {}

    template <typename>
    friend std::ostream& operator<<(std::ostream& out, const Disk& p);

    point_type center;
    double radius;
  };

  template <typename coord_type>
  std::ostream& operator<<(std::ostream& out, const Disk<coord_type>& p) {
    out << "[ center: (" << p.center << "), radius: " << p.radius << " ]";
    return out;
  }


  template <typename> class ConvexPolygon;

  template <typename CoordT>
  class Polygon : public std::vector<Point<CoordT> > {
  public:
    typedef Point<CoordT>       point_type;
    typedef Segment<CoordT>     segment_type;
    typedef BoundingBox<CoordT> bounding_box_type;

    Polygon() : _bbMemo(false) {}

    const BoundingBox<CoordT>& boundingBox() const {
      if (!_bbMemo) {
        CoordT xmin, xmax, ymin, ymax;
        typename Polygon<CoordT>::const_iterator point = this->begin();
        typename Polygon<CoordT>::const_iterator end = this->end();
        xmin = point->x(); xmax = point->x();
        ymin = point->y(); ymax = point->y();
        for ( ; point != end; ++point) {
          xmin = std::min(xmin, point->x());
          xmax = std::max(xmax, point->x());
          ymin = std::min(ymin, point->y());
          ymax = std::max(ymax, point->y());
        }
        _bbMemo = true;
        _bb = BoundingBox<CoordT>(xmin, ymin, xmax, ymax);
      }
      return _bb;
    }

    /// Graham's scan algorithm
    /// Cormen, Thomas H.; Leiserson, Charles E., Rivest, Ronald L., Stein, Clifford (2001) [1990].
    /// "33.3: Finding the convex hull".
    /// Introduction to Algorithms (2nd ed.).
    /// MIT Press and McGraw-Hill. pp. pp. 955–956. ISBN 0-262-03293-7.
  private:
    struct ccw {
      bool operator()(point_type i, point_type j) {
        return counterClockwiseTurn(p0, i, j) > 0;
      }
      point_type p0;
    };
  public:
    ConvexPolygon<CoordT> convexHull() const {
      // point with the lowest ordinate, or, in case of equality,
      // the lowest ordinate-leftest abscissa
      typename std::vector<point_type>::const_iterator it = this->begin();
      typename std::vector<point_type>::const_iterator lowest = it;
      for ( ; it != this->end(); ++it) {
        if (it->y() < lowest->y()
            || (it->y() == lowest->y() && it->x() < lowest->x()))
          lowest = it;
      }
      std::vector<point_type> points;
      std::copy(this->begin(), lowest, std::back_inserter(points));
      std::copy(lowest + 1, this->end(), std::back_inserter(points));
      // sort the points counter-clockwise relative to p0
      ccw comp;
      comp.p0 = *lowest;
      std::sort(points.begin(), points.end(), comp);
      // stack to the hull while turning left,
      // pop otherwise
      std::deque<point_type> hull;
      hull.push_back(*lowest);
      hull.push_back(points[0]);
      hull.push_back(points[1]);
      for (size_t i = 2; i < points.size(); ++i) {
        while (counterClockwiseTurn(hull[hull.size()-1],
                                    hull[hull.size()-2],
                                    points[i]) > 0)
          hull.pop_back();
        hull.push_back(points[i]); 
      }
      // return the polygon
      ConvexPolygon<CoordT> convexHull;
      for (typename std::deque<point_type>::iterator p = hull.begin(); p != hull.end(); ++p)
        convexHull.addPoint(*p);
      return convexHull;
    }

    void addPoint(const point_type& p) {
      this->push_back(p);

      _segments.clear();
      typename std::vector<point_type>::const_iterator it = this->begin();
      typename std::vector<point_type>::const_iterator first = it;
      typename std::vector<point_type>::const_iterator previous = it++;
      for (; it != this->end(); ++it) {
        _segments.push_back(segment_type(*previous, *it));
        previous = it;
      }
      _segments.push_back(segment_type(*previous, *first));
    }

    /// Iterations

    using typename std::vector<point_type>::const_iterator;

    typedef typename std::vector<segment_type>::const_iterator const_segment_iterator;
    const_segment_iterator segment_begin() const {
      return _segments.begin();
    }
    const_segment_iterator segment_end() const {
      return _segments.end();
    }

    template <typename>
    friend std::ostream& operator<<(std::ostream& out, const Polygon& p);

  private:
    // polygon as a closed set of vectors
    std::vector<segment_type> _segments;
    // memoize the bounding box
    mutable bool              _bbMemo;
    mutable bounding_box_type _bb;
  };

  template <typename CoordT>
  class ConvexPolygon : public Polygon<CoordT> {
  private:
    ConvexPolygon() {};
    template <typename> friend class Polygon;
  };

  template <typename coord_type>
  std::ostream& operator<<(std::ostream& out, const Polygon<coord_type>& p) {
    out << "[ " << p.size() << " vertices: (";
    typename Polygon<coord_type>::const_iterator vit;
    for (vit = p.begin(); vit != p.end(); ++vit)
      out << "(" << *vit << ')';
    out << ") ]";
    return out;
  }


  template <typename Region>
  class Union : public std::vector<Region> {
  public:
    typedef typename Region::point_type point_type;
    typedef typename std::vector<Region>::iterator iterator;
    typedef typename std::vector<Region>::const_iterator const_iterator;
    typedef Union<typename Region::bounding_box_type> bounding_box_type;
  };

  template <typename Region>
  class Intersection : public std::vector<Region> {
  public:
    typedef typename Region::point_type point_type;
    typedef typename std::vector<Region>::iterator iterator;
    typedef typename std::vector<Region>::const_iterator const_iterator;
  };

  template <typename Region1, typename Region2>
  class Difference : public std::pair<Region1,Region2> {
    typedef typename std::pair<Region1,Region2>::second_type complement;
  };

  template <typename Region>
  class Complement : public Region {
  };

  template <typename Region>
  struct has_center {
    static const bool value = false;
  };

  template <typename Coord>
  struct has_center<Disk<Coord> > {
    static const bool value = true;
  };


  /*************************************************
   *
   * Distances (Metrics).
   *
   *************************************************/
  namespace Euclid {
    /// norm of the cross-product
    template <typename Coord>
    inline
    Coord cross_product(Point<Coord> a, Point<Coord> b, Point<Coord> c) {
      return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
    }

    /// Euclidean distance, squared
    class distance : public std::binary_function<Point<int>,Point<int>,double> {
    public:
      template <typename Coord>
      double operator()(Point<Coord> from, Point<Coord> to) const {
        Coord x2 = from.x() - to.x();
        Coord y2 = from.y() - to.y();
        return std::sqrt(double(x2 * x2 + y2 * y2));
      }

      template <typename Coord>
      double operator()(Segment<Coord> segment, Point<Coord> point) const {
        double dist = Geo::Euclid::cross_product(segment.from(), segment.to(), point) / operator()(segment.from(), segment.to());
        Coord dot1 = dot_product(segment.from(), segment.to(), point);
        if (dot1 > 0)
          return operator()(segment.to(), point);
        Coord dot2 = dot_product(segment.to(), segment.from(), point);
        if (dot2 > 0)
          return operator()(segment.from(), point);
        return std::abs(dist);
      }
    };
  }


  namespace WGS84 {
    /// Convert the coordinates from fixed-precision numbers.
    ///
    /// Six digits are used for the floating part.
    /// This is enough for an accuracy of ~11cm at the equator.
    struct Point {
      Point(Geo::Point<int> p)
      : lat(double(p.x()) / 1000000),
        lon(double(p.y()) / 1000000) {}

      Point(double lat, double lon)
      : lat(lat), lon(lon) {}

      Point(Geo::Point<double> p)
      : lat(p.x()), lon(p.y()) {}

      double lat;
      double lon;
    };

    struct Segment {
      Segment(Geo::Segment<int> s)
      : from(s.from()),
        to(s.to()) {}

      Segment(Geo::Segment<double> s)
      : from(s.from()),
        to(s.to()) {}

      Point from;
      Point to;
    };

    /// Distance in meters for points representing WGS-84, the most used
    /// geodetic system.
    ///
    /// The distance function is Haversines, that models
    /// the earth as a sphere but is quite efficient.
    class distance : public std::binary_function<Point,Point,double> {
    public:
      double operator()(Point from, Point to) const {
        const double latRad = (from.lat - to.lat) * DEG_TO_RAD();
        const double lonRad = (from.lon - to.lon) * DEG_TO_RAD();
        const double latH = haversin(latRad);
        const double lonH = haversin(lonRad);

        const double tmp = cos(from.lat * DEG_TO_RAD()) * cos(to.lat * DEG_TO_RAD());
        const double arc = 2.0 * asin(std::sqrt(latH + tmp * lonH));
        return arc * EARTH_RADIUS_IN_METERS();
      }

    private:
      static double haversin(double angle) {
        double tmp = sin(angle * 0.5);
        return tmp * tmp;
      }

      static double DEG_TO_RAD() { return 0.01745329; }
      static double EARTH_RADIUS_IN_METERS() { return 6378137; }
    };

    /// If performance is an issue and accuracy less important,
    /// for small distances Pythagoras’ theorem can be used
    /// on an equirectangular projection.
    /// http://http://en.wikipedia.org/wiki/Equirectangular_projection
    class equirectangular_projection_distance
    : public std::binary_function<Point,Point,double> {
    public:
      double operator()(Point from, Point to) const {
        const double x = (from.lon - to.lon) * DEG_TO_RAD() * cos((from.lat + to.lat) / 2 * DEG_TO_RAD());
        const double y = (from.lat - to.lat) * DEG_TO_RAD();
        return std::sqrt(x * x + y * y) * EARTH_RADIUS_IN_METERS();
      }

      double operator()(Segment segment, Point point) const {
        const double cosphy = cos((segment.from.lat + segment.to.lat + point.lat) / 3 * DEG_TO_RAD());
        const double tmp = cosphy * DEG_TO_RAD();
        Geo::Segment<double> projected_segment(Geo::Point<double>(segment.from.lon * tmp, segment.from.lat * DEG_TO_RAD()),
                                               Geo::Point<double>(segment.to.lon * tmp, segment.to.lat * DEG_TO_RAD()));
        Geo::Point<double> projected_point(point.lon * tmp, point.lat * DEG_TO_RAD());
        return distance(projected_segment, projected_point) * EARTH_RADIUS_IN_METERS();
      }

    private:
      Geo::Euclid::distance distance;

      static double DEG_TO_RAD() { return 0.01745329; }
      static double EARTH_RADIUS_IN_METERS() { return 6378137; }
    };
  }


  /*************************************************
   *
   * Spatial predicates.
   *
   *************************************************/

  /// ---------------- Containment
  template <typename T1, typename T2, typename Distance>
  struct contains : public std::binary_function<T1,T2,bool> {
    bool operator()(const T1&, const T2&) const;
  };

  /// Disk specialization
  template <typename Coord, typename Distance>
  struct contains<Disk<Coord>,Point<Coord>,Distance>
  : public std::binary_function<Disk<Coord>,Point<Coord>,bool> {
    bool operator()(const Disk<Coord>& disk, const Point<Coord>& point) const {
      return disk.radius >= distance(disk.center, point);
    }
    Distance distance;
  };

  /// Polygon specialization
  /// Jordan curve theorem applied to polygons
  /// http://en.wikipedia.org/wiki/Jordan_curve_theorem
  /// Shoot a test ray along +X axis.
  /// The strategy is to compare vertex Y values to the testing point's Y and
  /// quickly discard edges which are entirely to one side of the test ray.
  /// From  Haines, Eric, "Point in Polygon Strategies," Graphics Gems IV,
  /// ed. Paul Heckbert, Academic Press, p. 24-46, 1994. 
  template <typename Coord, typename Distance>
  struct contains<Polygon<Coord>,Point<Coord>,Distance>
  : public std::binary_function<Polygon<Coord>,Point<Coord>,bool> {
    typedef Polygon<Coord> polygon_type;
    typedef Point<Coord>   point_type;

    bool operator()(const polygon_type& polygon, const point_type& point) const {
      // beware of overflows with integer types!
      double y = point.y();
      double x = point.x();
      bool inside_flag = false;
      // get test bool for above/below X axis
      typename polygon_type::const_segment_iterator seg = polygon.segment_begin();
      bool yflag0 = (seg->from().y() >= point.y());

      for ( ;seg != polygon.segment_end(); ++seg) {
        bool yflag1 = (seg->to().y() >= point.y());
        // Check if endpoints straddle (are on opposite sides) of X axis
        // (i.e. the Y's differ); if so, +X ray could intersect this edge.
        if (yflag0 != yflag1) {
          // Check intersection of segment with +X ray.
          // Note if >= point's X; if so, the ray hits it.
          // The division operation is avoided for the ">=" test by checking
          // the sign of the first vertex with respect to the test point.
          if (((seg->to().y() - y) * (seg->from().x() - seg->to().x()) >=
               (seg->to().x() - x) * (seg->from().y() - seg->to().y())) == yflag1) {
            inside_flag = !inside_flag;
          }
          yflag0 = yflag1;
        }
      }
      return inside_flag;
    }
  };

  /// Disk in Polygon specialization
  /// A disk is contained if its center is contained
  /// and it does not intersect with any edge.
  template <typename Coord, typename Distance>
  struct contains<Polygon<Coord>,Disk<Coord>,Distance>
  : public std::binary_function<Polygon<Coord>,Disk<Coord>,bool> {
    typedef Polygon<Coord> polygon_type;
    typedef Disk<Coord>    disk_type;
    typedef Point<Coord>   point_type;
    typedef Distance       distance_type;

    bool operator()(const polygon_type& polygon, const disk_type& disk) const {
      if (!_contains(polygon, disk.center))
        return false;

      // Determine the point from the polygon boundary that is the
      // closest to the ball's center
      double closest_distance_to_ball_center = std::numeric_limits<double>::max();
      for (typename polygon_type::const_segment_iterator seg = polygon.segment_begin();
          seg != polygon.segment_end(); ++seg) {
        closest_distance_to_ball_center = std::min(closest_distance_to_ball_center,
            distance(*seg, disk.center));
      }
      return closest_distance_to_ball_center >= disk.radius;
    }

  private:
    contains<polygon_type,point_type,distance_type> _contains;
    distance_type distance;
  };

  // force the use of the equirectangular projection in WGS84
  template <typename Coord>
  struct contains<Polygon<Coord>,Disk<Coord>,Geo::WGS84::distance>
  : public std::binary_function<Polygon<Coord>,Disk<Coord>,bool> {
    bool operator()(const Polygon<Coord>& polygon, const Disk<Coord>& d) const {
      return f(polygon, d);
    }

    contains<Polygon<Coord>,
      Disk<Coord>,
      Geo::WGS84::equirectangular_projection_distance> f;
  };

  /// Bounding Box specialization
  template <typename Coord, typename Distance>
  struct contains<BoundingBox<Coord>,Point<Coord>,Distance>
  : public std::binary_function<BoundingBox<Coord>,Point<Coord>,bool> {
    bool operator()(const BoundingBox<Coord>& bb, const Point<Coord>& point) const {
      return bb.xMin() <= point.x() && point.x() <= bb.xMax() &&
        bb.yMin() <= point.y() && point.y() <= bb.yMax();
    }
  };

  /// Bounding Box inside Bounding Box
  template <typename Coord, typename Distance>
  struct contains<BoundingBox<Coord>,BoundingBox<Coord>,Distance>
  : public std::binary_function<BoundingBox<Coord>,BoundingBox<Coord>,bool> {
    bool operator()(const BoundingBox<Coord>& bb1, const BoundingBox<Coord>& bb2) const {
      return bb1.xMin() <= bb2.xMin() && bb1.yMin() <= bb2.yMin()
        && bb1.xMax() >= bb2.xMax() && bb1.yMax() >= bb2.yMax();
    }
  };

  /// Union specialization
  template <typename Region1, typename Region2, typename Distance>
  struct contains<Union<Region1>,Region2,Distance>
  : public std::binary_function<Union<Region1>,Region2,bool> {
    bool operator()(const Union<Region1>& uni, const Region2& thing) const {
      for (typename Union<Region1>::const_iterator r = uni.begin();
          r != uni.end(); ++r) {
        if (contains<Region1,Region2,Distance>()(*r, thing))
            return true;
      }
      return false;
    }
  };

  /// Intersection specialization
  template <typename Region1, typename Region2, typename Distance>
  struct contains<Intersection<Region1>,Region2,Distance>
  : public std::binary_function<Intersection<Region1>,Region2,bool> {
    bool operator()(const Intersection<Region1>& inter, const Region2& thing) const {
      if (inter.empty()) return false;
      for (typename Intersection<Region1>::const_iterator r = inter.begin();
          r != inter.end(); ++r) {
        if (!contains<Region1,Region2,Distance>()(*r, thing))
          return false;
      }
      return true;
    }
  };

  /// Difference specialization
  template <typename Region1, typename Region2, typename Region3, typename Distance>
  struct contains<Difference<Region1,Region2>,Region3,Distance>
  : public std::binary_function<Difference<Region1,Region2>,Region3,bool> {
    bool operator()(const Difference<Region1,Region2>& diff, const Region3& bb) const {
      return contains<Region1,Region3,Distance>()(diff.first, bb)
        && ! contains<Region2,Region3,Distance>()(diff.second, bb);
    }
  };

  /// Complement specialization
  template <typename Region1, typename Region2, typename Distance>
  struct contains<Complement<Region1>,Region2,Distance>
  : public std::binary_function<Complement<Region1>,Region2,bool> {
    bool operator()(const Complement<Region1>& complem, const Region2& bb) const {
      return !contains<Region1,Region2,Distance>()(complem, bb);
    }
  };



  /// ---------------- Intersection
  template <typename T1, typename T2, typename Distance>
  struct intersects
  : public commutative_binary_function<intersects<T1,T2,Distance>,T1,T2,bool>
  {};

  /// Bounding Box - Bounding Box specialization
  /// Hyperplane separation theorem
  /// http://http://en.wikipedia.org/wiki/Hyperplane_separation_theorem
  template <typename Coord, typename Distance>
  struct intersects<BoundingBox<Coord>,BoundingBox<Coord>,Distance>
  : public std::binary_function<BoundingBox<Coord>,BoundingBox<Coord>,bool> {
    typedef BoundingBox<Coord> bounding_box_type;

    bool operator()(const bounding_box_type& bb1, const bounding_box_type& bb2) const {
      return !(bb1.xMin() > bb2.xMax()
          || bb1.xMax() < bb2.xMin()
          || bb1.yMax() < bb2.yMin()
          || bb1.yMin() > bb2.yMax());
    }
  };

  /// Convex Polygon - Bounding Box specialization
  /// Hyperplane separation theorem
  /// http://http://en.wikipedia.org/wiki/Hyperplane_separation_theorem
  template <typename Coord, typename Distance>
  struct intersects<ConvexPolygon<Coord>,BoundingBox<Coord>,Distance>
  : public std::binary_function<ConvexPolygon<Coord>,BoundingBox<Coord>,bool> {
    typedef BoundingBox<Coord>    bounding_box_type;
    typedef ConvexPolygon<Coord>  polygon_type;
    typedef Point<Coord>          point_type;
    typedef Coord                 coord_type;

    bool operator()(const polygon_type& polygon, const bounding_box_type& bb) const {
      const bounding_box_type& pbb = polygon.boundingBox();
      if (!_intersects(pbb, bb))
        return false; // separating line along an axis

      // Project the bounding box along the perpendicular axis to each
      // polygon's edge.
      for (typename polygon_type::const_segment_iterator seg = polygon.segment_begin();
          seg != polygon.segment_end(); ++seg) {
        // perpendicular axis passing through a endpoint of the segment
        point_type v(seg->from().x() + seg->from().y() - seg->to().y(),
                     seg->from().y() + seg->to().x() - seg->from().x());
        // projection of the bounding box on the axis
        coord_type minBB = dot_product(seg->from(), v, bb.first);
        coord_type maxBB = minBB;
        coord_type dot2 = dot_product(seg->from(), v, bb.second);
        minBB = std::min(minBB, dot2);
        maxBB = std::max(maxBB, dot2);
        coord_type dot3 = dot_product(seg->from(), v, point_type(bb.xMin(), bb.yMax()));
        minBB = std::min(minBB, dot3);
        maxBB = std::max(maxBB, dot3);
        coord_type dot4 = dot_product(seg->from(), v, point_type(bb.xMax(), bb.yMin()));
        minBB = std::min(minBB, dot4);
        maxBB = std::max(maxBB, dot4);
        // projection of the polygon on the axis
        typename polygon_type::const_iterator point = polygon.begin();
        coord_type minP = dot_product(seg->from(), v, *point);
        coord_type maxP = minP;
        for ( ; point != polygon.end(); ++point) {
          coord_type dp = dot_product(seg->from(), v, *point);
          minP = std::min(minP, dp);
          maxP = std::max(maxP, dp);
        }
        // are the projections overlapping?
        if (minBB > maxP || maxBB < minP)
          return false;
      }
      return true;
    }

  private:
    intersects<bounding_box_type,bounding_box_type,Distance> _intersects;
  };

  /// Disk - Disk specialization
  template <typename Coord, typename Distance>
  struct intersects<Disk<Coord>,Disk<Coord>,Distance>
  : public std::binary_function<Disk<Coord>,Disk<Coord>,bool> {
    typedef Disk<Coord>        ball_type;
    typedef Point<Coord>       point_type;
    typedef Distance           distance_type;

    bool operator()(const ball_type& b1, const ball_type& b2) const {
      return distance(b1.center, b2.center) <= b1.radius + b2.radius;
    }

    distance_type distance;
  };

  /// Disk - Bounding Box specialization
  template <typename Coord, typename Distance>
  struct intersects<Disk<Coord>,BoundingBox<Coord>,Distance>
  : public std::binary_function<Disk<Coord>,BoundingBox<Coord>,bool> {
    typedef BoundingBox<Coord> bounding_box_type;
    typedef Disk<Coord>        ball_type;
    typedef Point<Coord>       point_type;
    typedef Coord              coord_type;
    typedef Distance           distance_type;

    bool operator()(const ball_type& d, const bounding_box_type bb) const {
      if (_contains(bb, d.center)) {
        return true;
      }

      // Determine the point from the bounding box boundary that is the
      // closest to the ball's center (9 possible cases)
      coord_type closestX;
      if (bb.xMin() <= d.center.x() && d.center.x() <= bb.xMax()) {
        closestX = d.center.x();
      } else {
        if (std::abs(bb.xMin() - d.center.x()) > std::abs(bb.xMax() - d.center.x())) {
          closestX = bb.xMax();
        } else {
          closestX = bb.xMin();
        }
      }

      coord_type closestY;
      if (bb.yMin() <= d.center.y() && d.center.y() <= bb.yMax()) {
        closestY = d.center.y();
      } else {
        if (std::abs(bb.yMin() - d.center.y()) > std::abs(bb.yMax() - d.center.y())) {
          closestY = bb.yMax();
        } else {
          closestY = bb.yMin();
        }
      }

      point_type closestFromCenter(closestX, closestY);
      return distance(d.center, closestFromCenter) <= d.radius;
    }

    distance_type distance;
    contains<bounding_box_type,point_type,distance_type> _contains;
  };

  /// Polygon - Disk specialization
  /// In spherical geometry, use the equirectangular projection
  /// to be able to employ the Pythagoras theorem.
  template <typename Coord, typename Distance>
  struct intersects<Polygon<Coord>,Disk<Coord>,Distance>
  : public std::binary_function<Polygon<Coord>,Disk<Coord>,bool> {
    typedef Polygon<Coord>     polygon_type;
    typedef Disk<Coord>        ball_type;
    typedef Point<Coord>       point_type;
    typedef Distance         distance_type;

    bool operator()(const polygon_type& polygon, const ball_type& d) const {
      if (contains<polygon_type,point_type,distance_type>()(polygon, d.center)) {
        return true;
      }

      // Determine the point from the polygon boundary that is the
      // closest to the ball's center
      double closest_distance_to_ball_center = std::numeric_limits<double>::max();
      for (typename polygon_type::const_segment_iterator seg = polygon.segment_begin();
          seg != polygon.segment_end(); ++seg) {
        closest_distance_to_ball_center = std::min(closest_distance_to_ball_center,
                                                   distance(*seg, d.center));
      }
      return closest_distance_to_ball_center <= d.radius;
    }

    distance_type distance;
  };

  // force the use of the equirectangular projection in WGS84
  template <typename Coord>
  struct intersects<Polygon<Coord>,Disk<Coord>,Geo::WGS84::distance>
  : public std::binary_function<Polygon<Coord>,Disk<Coord>,bool> {
    bool operator()(const Polygon<Coord>& polygon, const Disk<Coord>& d) const {
      return f(polygon, d);
    }

    intersects<Polygon<Coord>,
               Disk<Coord>,
               Geo::WGS84::equirectangular_projection_distance> f;
  };

  /// Union specialization
  template <typename Region1, typename Region2, typename Distance>
  struct intersects<Union<Region1>,Region2,Distance>
  : public std::binary_function<Union<Region1>,Region2,bool> {
    bool operator()(const Union<Region1>& uni, const Region2& bb) const {
      for (typename Union<Region1>::const_iterator r = uni.begin();
          r != uni.end(); ++r) {
        if (intersects<Region1,Region2,Distance>()(*r, bb))
            return true;
      }
      return false;
    }
  };

  /// Intersection specialization
  template <typename Region1, typename Region2, typename Distance>
  struct intersects<Intersection<Region1>,Region2,Distance>
  : public std::binary_function<Intersection<Region1>,Region2,bool> {
    bool operator()(const Intersection<Region1>& inter, const Region2& bb) const {
      if (inter.empty()) return false;
      for (typename Intersection<Region1>::const_iterator r = inter.begin();
          r != inter.end(); ++r) {
        if (!intersects<Region1,Region2,Distance>()(*r, bb))
          return false;
      }
      return true;
    }
  };

  /// Difference specialization
  template <typename Region1, typename Region2, typename Region3, typename Distance>
  struct intersects<Difference<Region1,Region2>,Region3,Distance>
  : public std::binary_function<Difference<Region1,Region2>,Region3,bool> {
    bool operator()(const Difference<Region1,Region2>& diff, const Region3& bb) const {
      return intersects<Region1,Region3,Distance>()(diff.first, bb)
        && ! contains<Region2,Region3,Distance>()(diff.second, bb);
    }
  };

  /// Complement specialization
  template <typename Region1, typename Region2, typename Distance>
  struct intersects<Complement<Region1>,Region2,Distance>
  : public std::binary_function<Complement<Region1>,Region2,bool> {
    bool operator()(const Complement<Region1>& complem, const Region2& bb) const {
      return !intersects<Region1,Region2,Distance>()(complem, bb);
    }
  };

} // namespace Geometry

#endif
