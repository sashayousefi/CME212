#pragma once
#ifndef CME212_SFML_VIEWER_HPP
#define CME212_SFML_VIEWER_HPP

/**
 * @file SFML_Viewer.hpp
 * Interactive OpenGL Viewer
 *
 * @brief Provides a OpenGL window with some basic interactivity.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <map>
#include <numeric>

// Suppress compiler warnings generated from OpenGL deprecation on Mac OS
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <SFML/Window.hpp>
#include <SFML/OpenGL.hpp>

#include "GLCamera.hpp"
#include "Color.hpp"
#include "Point.hpp"
#include "Font5x7.hpp"

#ifdef __linux__
#include <X11/Xlib.h>
struct XInit { XInit() { XInitThreads(); } };
#endif

namespace CME212 {

// Map C++ types to GL types at compilation time: gltype<T>::value
template <GLenum V> struct gltype_v { static constexpr GLenum value = V; };
template <class T> struct gltype {};
template<> struct gltype<unsigned char>  : gltype_v<GL_UNSIGNED_BYTE>  {};
template<> struct gltype<unsigned short> : gltype_v<GL_UNSIGNED_SHORT> {};
template<> struct gltype<unsigned int>   : gltype_v<GL_UNSIGNED_INT>   {};
template<> struct gltype<short>          : gltype_v<GL_SHORT>          {};
template<> struct gltype<int>            : gltype_v<GL_INT>            {};
template<> struct gltype<float>          : gltype_v<GL_FLOAT>          {};
template<> struct gltype<double>         : gltype_v<GL_DOUBLE>         {};

// A default color functor that returns white for anything it recieves
struct DefaultColor {
  template <typename NODE>
  Color operator()(const NODE&) {
    return Color(1);
  }
};

// A default position functor that returns node.position() for any node
struct DefaultPosition {
  template <typename NODE>
  Point operator()(const NODE& node) {
    return node.position();
  }
};

/** Print any outstanding OpenGL error to std::cerr. */
void check_gl_error(const char* context = nullptr) {
  GLenum errCode = glGetError();;
  if (errCode != GL_NO_ERROR) {
    std::cerr << "OpenGL Error";
    if (context)
      std::cerr << " at " << context;
    std::cerr << ": " << errCode << "\n";
  }
}

/** Class that provides a graphical window context for drawing points and edges.
 */
class SFML_Viewer
{
  using safe_lock = std::lock_guard<std::mutex>;

 public:

  /** Constructor */
  SFML_Viewer(int sx = 640, int sy = 480)
      : window_(sf::VideoMode(sx, sy), "CME212",
                sf::Style::Resize | sf::Style::Close,
                sf::ContextSettings(24 /*depth*/, 8 /*stencil*/, 2 /*anti*/)),
        render_requested_(true),
        mouse_mask_(0)
  {
    init();
    render_thread_ = std::thread([&]{ render_loop(); });
  }

  /** Disable copy and assignment of SFML_Viewers */
  SFML_Viewer(const SFML_Viewer&) = delete;
  void operator=(const SFML_Viewer&) = delete;

  ~SFML_Viewer()
  {
    render_flag_.notify_one();
    render_thread_.join();
  }

  /** Main Event Loop
   *
   * Executed by the event thread until interrupt or killed
   */
  void event_loop()
  {
    sf::Event event;
    while (window_.waitEvent(event)) {
      safe_lock lock(mutex_);
      handle_event(event);
    }
  }

  /** Add the node positions in @a g to the graphics.
   * @pre The G type and @a g object must have the following properties:
   *   <ol><li>G::size_type and G::node_type are types.</li>
   *       <li>G::size_type g.num_nodes() returns the number of nodes.</li>
   *	     <li>G::node_type g.node(G::size_type i) returns a node.</li>
   *       <li>Point g.node(i).position() returns the position of node i.</li>
   *   </ol>
   *
   * This draws only nodes. All nodes are drawn as white.
   */
  template <typename G>
  void draw_graph_nodes(const G& g)
  {
    // Lock for data update
    { safe_lock lock(mutex_);

      // Convenience aliases
      using size_type = typename G::size_type;

      // Insert the nodes
      size_type num_nodes = g.num_nodes();
      for (size_type i = 0; i < num_nodes; ++i)
        coords_.push_back(g.node(i).position());

      // Resize with all added nodes white
      colors_.resize(coords_.size(), Color(1,1,1));
    }

    request_render();
  }

  /** Add the node positions and their edges in @a g to the graphics.
   * @pre The G type and @a g object must have the properties listed in
   *   draw_graph_nodes(), and in addition:
   *   <ol><li>G::edge_type is a type.</li>
   *       <li>G::size_type g.node(i).index() returns the node index.</li>
   *       <li>G::size_type g.num_edges() returns the number of edges.</li>
   *       <li>G::edge_type g.edge(G::size_type i) returns an edge.</li>
   *       <li>G::node_type g.edge(i).node1() and node2() return nodes.</li>
   *   </ol>
   *
   * This draws both nodes and edges. All nodes/edges are drawn as white.
   */
  template <typename G>
  void draw_graph(const G& g)
  {
    // Lock for data update
    { safe_lock lock(mutex_);

      // Convenience aliases
      using size_type = typename G::size_type;
      using node_type = typename G::node_type;
      using edge_type = typename G::edge_type;

      // Map needed if coords_ and/or edges_ already has data
      std::map<size_type, unsigned> nodemap;

      // Insert the nodes and record mapping
      size_type num_nodes = g.num_nodes();
      for (size_type i = 0; i < num_nodes; ++i) {
        node_type n = g.node(i);
        nodemap[n.index()] = coords_.size();
        coords_.push_back(n.position());
      }

      // Resize with all added nodes white
      colors_.resize(coords_.size(), Color(1,1,1));

      // Insert edges according to our mapping
      size_type num_edges = g.num_edges();
      for (size_type i = 0; i < num_edges; ++i) {
        edge_type e = g.edge(i);
        auto it1 = nodemap.find(e.node1().index());
        auto it2 = nodemap.find(e.node2().index());
        if (it1 != nodemap.end() && it2 != nodemap.end()) {
          edges_.push_back(it1->second);
          edges_.push_back(it2->second);
        }
      }
    }

    request_render();
  }

  /** Return an empty node map for the input graph.
   * @param  Not used. Convenience for type deduction.
   * @tparam G Type that defines a node_type:  typename G::node_type
   *
   * Node maps are passed to, and modified by, add_nodes() and add_edges().
   */
  template <typename G>
  std::map<typename G::node_type, unsigned> empty_node_map(const G&) const
  {
    return std::map<typename G::node_type, unsigned>();
  }

  /** Add the nodes in the range [first, last) to the display.
   * @param[in] first,last   Iterator range of nodes
   * @param[in,out] node_map Tracks node identities for use by add_edges().
   *    Create this argument by calling empty_node_map(graph).
   * @tparam InputIterator   Iterator type over nodes.
   *     std::is_same<typename std::iterator_traits<InputIterator>::value_type,
   *                  typename Map::key_type>::value
   *     requires value_type (node) to provide   Point ::position() const
   * @tparam Map     std::map with the key_type satisfying the above.
   *
   * @post For all i in [@a first, @a last), node_map.count(*i) == 1.
   *
   * Uses white for color and node.position() for position. It's OK to add a
   * Node more than once. Second and subsequent adds update the existing
   * node's position.
   */
  template <typename InputIterator, typename Map>
  void add_nodes(InputIterator first, InputIterator last,
                 Map& node_map)
  {
    return add_nodes(first, last, DefaultColor(), DefaultPosition(), node_map);
  }

  /** Add the nodes in the range [first, last) to the display.
   * @param[in] first,last    Iterator range of nodes
   * @param[in] color_functon Function from a node to a color.
   *                            CME212::Color color_function(node)
   * @param[in,out] node_map Tracks node identities for use by add_edges().
   *    Create this argument by calling empty_node_map(graph).
   * @tparam InputIterator   Iterator type over nodes.
   *     std::is_same<typename std::iterator_traits<InputIterator>::value_type,
   *                  typename Map::key_type>::value
   *     requires value_type (node) to provide   Point ::position() const
   * @tparam ColorFn Function type that maps value_type (node) to CME212:Color.
   * @tparam Map     std::map with the key_type satisfying the above.
   *
   * @post For all i in [@a first, @a last), node_map.count(*i) == 1.
   *
   * Uses color_function for color and node.position() for position. It's OK
   * to add a Node more than once. Second and subsequent adds update the
   * existing node's position.
   */
  template <typename InputIterator, typename ColorFn, typename Map>
  void add_nodes(InputIterator first, InputIterator last,
                 ColorFn color_function, Map& node_map)
  {
    return add_nodes(first, last, color_function, DefaultPosition(), node_map);
  }

  /** Add the nodes in the range [first, last) to the display.
   * @param[in] first,last    Iterator range of nodes
   * @param[in] color_functon Function from a node to a color.
   *                            CME212::Color color_function(node)
   * @param[in] position_function Function from a node to a position.
   *                            Point position_function(node)
   * @param[in,out] node_map Tracks node identities for use by add_edges().
   *    Create this argument by calling empty_node_map(graph).
   * @tparam InputIterator   Iterator type over nodes.
   *     std::is_same<typename std::iterator_traits<InputIterator>::value_type,
   *                  typename Map::key_type>::value
   * @tparam ColorFn Function type that maps value_type (node) to CME212:Color.
   * @tparam PointFn Function type that maps value_type (node) to Point.
   * @tparam Map     std::map with the key_type satisfying the above.
   *
   * @post For all i in [@a first, @a last), node_map.count(*i) == 1.
   *
   * Uses color_function for color and position_function for position. It's OK
   * to add a Node more than once. Second and subsequent adds update the
   * existing node's position.
   */
  template <typename InputIterator,
            typename ColorFn, typename PointFn, typename Map>
  void add_nodes(InputIterator first, InputIterator last,
                 ColorFn color_function, PointFn position_function,
                 Map& node_map)
  {
    // Lock for data update
    { safe_lock lock(mutex_);

      for ( ; first != last; ++first) {
        // Get node and record the index mapping
        auto n = *first;
        auto r = node_map.insert(typename Map::value_type(n,coords_.size()));
        if (r.second) {   // new node was inserted
          coords_.push_back(position_function(n));
          colors_.push_back(color_function(n));
        } else {          // node already exists and not updated
          unsigned index = r.first->second;
          coords_[index] = position_function(n);
          colors_[index] = color_function(n);
        }
      }
    }

    request_render();
  }

  /** Add the edges in the range [first, last) to the display.
   * @param[in] first,last Iterator range of edges
   * @param[in] node_map Tracks node identities.
   *
   * @tparam InputIterator Iterator type over edges. An edge is a type
   *     that provides:
   *        using node_type = ...
   *        node_type node1()
   *        node_type node2()
   *     where std::is_same<node_type, typename Map::key_type>::value.
   *
   * Edges whose endpoints weren't previously added to the node_map by
   * add_nodes() (i.e. node_map.count(node) == 0) are ignored.
   */
  template <typename InputIterator, typename Map>
  void add_edges(InputIterator first, InputIterator last, const Map& node_map)
  {
    // Lock for data update
    { safe_lock lock(mutex_);

      for ( ; first != last; ++first) {
        auto edge = *first;
        auto n1 = node_map.find(edge.node1());
        auto n2 = node_map.find(edge.node2());
        if (n1 != node_map.end() && n2 != node_map.end()) {
          edges_.push_back(n1->second);
          edges_.push_back(n2->second);
        }
      }
    }

    request_render();
  }

  /** Set a string label to display "green LCD" style. */
  void set_label(const std::string& str)
  {
    safe_lock lock(mutex_);
    if (str != label_) {
      label_ = str;
      request_render();
    }
  }

  /** Set a label to display "green LCD" style. */
  void set_label(double d)
  {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4) << d;
    set_label(ss.str());
  }

  /** Center view.
   *
   * Attempts to center the OpenGL view on the object by setting the new
   * viewpoint to the average of all nodes
   */
  void center_view()
  {
    { safe_lock lock(mutex_);
      center_camera();
    }
    request_render();
  }

  /** Request that the screen update shortly. */
  void request_render()
  {
    if (!render_requested_) {
      render_requested_ = true;
      render_flag_.notify_one();
    }
  }

  /** Erase graphics. */
  void clear()
  {
    { safe_lock lock(mutex_);
      coords_.clear();
      colors_.clear();
      edges_.clear();
      label_.clear();
    }
    request_render();
  }

 private:

  /** Initialize the Window state
   */
  void init()
  {
    window_.setActive(true);

    // Background Color
    glClearColor(0.0, 0.0, 0.0, 1.0);
    // Initialize View
    glViewport(0, 0, window_.getSize().x, window_.getSize().y);

    // Set projection matrix
    camera_.set_perspective(60, window_.getSize().x/double(window_.getSize().y),
                            0.05, 1000);

    // Set up the camera view
    camera_.zoom_mag(2);
    //camera_.rotate_y(0.5);
    //camera_.rotate_x(0.5);

    // Point system
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // OpenGL Fog for better depth perception
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, GL_EXP);
    glFogf(GL_FOG_DENSITY, 0.25);

    window_.setActive(false);
  }

  /** Main Event Handler
   *
   * @param[in] event The mouse, keyboard, or screen event to be handled
   */
  void handle_event(sf::Event& event)
  {
    switch(event.type)
    {
      case sf::Event::MouseButtonPressed: {
        if (event.mouseButton.button == sf::Mouse::Right)
          mouse_mask_ |= BUTTON_RIGHT;
        if (event.mouseButton.button == sf::Mouse::Left)
          mouse_mask_ |= BUTTON_LEFT;
        last_mouse_x_ = event.mouseButton.x;
        last_mouse_y_ = event.mouseButton.y;
      } break;

      case sf::Event::MouseButtonReleased: {
        if (event.mouseButton.button == sf::Mouse::Right)
          mouse_mask_ &= ~BUTTON_RIGHT;
        if (event.mouseButton.button == sf::Mouse::Left)
          mouse_mask_ &= ~BUTTON_LEFT;
        last_mouse_x_ = event.mouseButton.x;
        last_mouse_y_ = event.mouseButton.y;
      } break;

      // The mouse moved over the screen
      case sf::Event::MouseMoved: {
        int xrel = event.mouseMove.x - last_mouse_x_;
        int yrel = event.mouseMove.y - last_mouse_y_;
        // Left mouse button is down
        if (mouse_mask_ & BUTTON_LEFT) {
          camera_.rotate_x( 0.01*yrel);
          camera_.rotate_y(-0.01*xrel);
          request_render();
        }
        // Right mouse button is down
        if (mouse_mask_ & BUTTON_RIGHT) {
          camera_.pan(-0.004*xrel, 0.004*yrel, 0);
          request_render();
        }
        last_mouse_x_ = event.mouseMove.x;
        last_mouse_y_ = event.mouseMove.y;
        // Avoid rendering on every mouse motion event
      } break;

      case sf::Event::MouseWheelScrolled: {
        // Wheel event
        if (event.mouseWheelScroll.delta > 0)
          camera_.zoom(1.10);
        if (event.mouseWheelScroll.delta < 0)
          camera_.zoom(0.91);

        request_render();
      } break;

      case sf::Event::KeyPressed: {
        // Keyboard 'c' to center
        if (event.key.code == sf::Keyboard::C) {
          center_camera();
          request_render();
        }
        // Keyboard 'esc' or 'q' to exit
        if (event.key.code == sf::Keyboard::Escape ||
            event.key.code == sf::Keyboard::Q)
          window_.close();
      } break;

      case sf::Event::Resized: {
        //resize();   // XXX: TODO
        request_render();
      } break;

      case sf::Event::Closed: {
        window_.close();
      } break;

      default:
        break;
    }
  }

  void center_camera()
  {
    // Sum the points in coords_
    Point c = std::accumulate(coords_.begin(), coords_.end(), Point(0));
    // Set the new view point
    camera_.view_point(c / coords_.size());
  }

  /** Render a label in "Green LCD" style. Only knows digits, spaces, and '.'
   */
  void render_label()
  {
    static const GLfloat skewscalem[] = {
      3, 0, 0, 0, 0.5, 3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

    // Set both relevant matrices for 2D display.
    // The projection matrix is orthographic.
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, window_.getSize().x, 0, window_.getSize().y, -1, 1);

    // The model view matrix is the identity.
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Draw in green with fat points
    glColor3f(0, .759, .437);
    glPointSize(3);

    double expected_width = label_.size()*5*3 + 15;
    // Translate and skew.
    glTranslatef(window_.getSize().x - expected_width, 10, 0);
    glMultMatrixf(skewscalem);

    // Draw label
    for (char label_c : label_) {
      const char* code = Font5x7 + (label_c - 32)*5;
      glBegin(GL_POINTS);
      for (char x = 0; x < 5; ++x)
        for (char c = *code++, y = 6; c; c >>= 1, --y)
          if (c & 1) glVertex2f(x, y);
      glEnd();
      glTranslatef(5, 0, 0);
    }

    // Restore settings.
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    check_gl_error();
  }

  /** Render to the screen
   *
   * Double buffered rendering using OpenGL
   */
  void render_loop()
  {
    std::unique_lock<std::mutex> lock(mutex_);

    do {
      // Activate the window for OpenGL rendering
      window_.setActive(true);

      // Clear the screen and z-buffer
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // Construct the view matrix
      camera_.set_view();

      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_COLOR_ARRAY);

      // Define the vertex interpreter
      glVertexPointer(Point::size(), gltype<typename Point::value_type>::value,
                      0, coords_.data());

      // Define the color interpreter
      glColorPointer(3, gltype<typename Color::value_type>::value,
                     0, colors_.data());

      // Draw the points
      glPointSize(1.5);
      glDrawArrays(GL_POINTS, 0, coords_.size());

      // Draw the lines
      glLineWidth(1);
      glDrawElements(GL_LINES, edges_.size(), gltype<unsigned>::value,
                     edges_.data());

      glDisableClientState(GL_COLOR_ARRAY);
      glDisableClientState(GL_VERTEX_ARRAY);

      // Draw the label
      if (!label_.empty())
        render_label();

      // Make visible
      window_.display();

      // Deactivate the window for OpenGL rendering
      window_.setActive(false);

      // Wait for next render request
      render_requested_ = false;
      render_flag_.wait(lock);

    } while (window_.isOpen());
  }

#ifdef __linux__
  XInit init_;   // Must be first member
#endif
  sf::Window window_;

  std::mutex mutex_;
  std::thread render_thread_;
  std::condition_variable render_flag_;
  bool render_requested_;

  static constexpr char BUTTON_RIGHT = (1 << 0);
  static constexpr char BUTTON_LEFT  = (1 << 1);
  char mouse_mask_;
  int last_mouse_x_, last_mouse_y_;

  // OpenGL Camera to track the current view
  GLCamera camera_;

  // Vertices
  std::vector<Point> coords_;

  // Colors
  std::vector<Color> colors_;

  // Edges (index pairs)
  std::vector<unsigned> edges_;

  // Currently displayed label
  std::string label_;
};

} // end namespace CME212

#endif
