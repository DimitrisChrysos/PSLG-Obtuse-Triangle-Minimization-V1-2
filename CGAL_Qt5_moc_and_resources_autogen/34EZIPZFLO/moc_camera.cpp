/****************************************************************************
** Meta object code from reading C++ file 'camera.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../../../../../usr/local/include/CGAL/Qt/camera.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'camera.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_CGAL__qglviewer__Camera_t {
    QByteArrayData data[75];
    char stringdata0[866];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_CGAL__qglviewer__Camera_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_CGAL__qglviewer__Camera_t qt_meta_stringdata_CGAL__qglviewer__Camera = {
    {
QT_MOC_LITERAL(0, 0, 23), // "CGAL::qglviewer::Camera"
QT_MOC_LITERAL(1, 24, 11), // "setPosition"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 3), // "Vec"
QT_MOC_LITERAL(4, 41, 3), // "pos"
QT_MOC_LITERAL(5, 45, 14), // "setOrientation"
QT_MOC_LITERAL(6, 60, 10), // "Quaternion"
QT_MOC_LITERAL(7, 71, 1), // "q"
QT_MOC_LITERAL(8, 73, 5), // "theta"
QT_MOC_LITERAL(9, 79, 3), // "phi"
QT_MOC_LITERAL(10, 83, 11), // "setUpVector"
QT_MOC_LITERAL(11, 95, 2), // "up"
QT_MOC_LITERAL(12, 98, 6), // "noMove"
QT_MOC_LITERAL(13, 105, 16), // "setViewDirection"
QT_MOC_LITERAL(14, 122, 9), // "direction"
QT_MOC_LITERAL(15, 132, 6), // "lookAt"
QT_MOC_LITERAL(16, 139, 6), // "target"
QT_MOC_LITERAL(17, 146, 15), // "showEntireScene"
QT_MOC_LITERAL(18, 162, 9), // "fitSphere"
QT_MOC_LITERAL(19, 172, 6), // "center"
QT_MOC_LITERAL(20, 179, 6), // "radius"
QT_MOC_LITERAL(21, 186, 14), // "fitBoundingBox"
QT_MOC_LITERAL(22, 201, 3), // "min"
QT_MOC_LITERAL(23, 205, 3), // "max"
QT_MOC_LITERAL(24, 209, 15), // "fitScreenRegion"
QT_MOC_LITERAL(25, 225, 9), // "rectangle"
QT_MOC_LITERAL(26, 235, 11), // "centerScene"
QT_MOC_LITERAL(27, 247, 24), // "interpolateToZoomOnPixel"
QT_MOC_LITERAL(28, 272, 5), // "pixel"
QT_MOC_LITERAL(29, 278, 21), // "interpolateToFitScene"
QT_MOC_LITERAL(30, 300, 13), // "interpolateTo"
QT_MOC_LITERAL(31, 314, 5), // "Frame"
QT_MOC_LITERAL(32, 320, 2), // "fr"
QT_MOC_LITERAL(33, 323, 8), // "duration"
QT_MOC_LITERAL(34, 332, 7), // "setType"
QT_MOC_LITERAL(35, 340, 4), // "Type"
QT_MOC_LITERAL(36, 345, 4), // "type"
QT_MOC_LITERAL(37, 350, 14), // "setFieldOfView"
QT_MOC_LITERAL(38, 365, 3), // "fov"
QT_MOC_LITERAL(39, 369, 24), // "setHorizontalFieldOfView"
QT_MOC_LITERAL(40, 394, 4), // "hfov"
QT_MOC_LITERAL(41, 399, 16), // "setFOVToFitScene"
QT_MOC_LITERAL(42, 416, 14), // "setAspectRatio"
QT_MOC_LITERAL(43, 431, 6), // "aspect"
QT_MOC_LITERAL(44, 438, 23), // "setScreenWidthAndHeight"
QT_MOC_LITERAL(45, 462, 5), // "width"
QT_MOC_LITERAL(46, 468, 6), // "height"
QT_MOC_LITERAL(47, 475, 16), // "devicePixelRatio"
QT_MOC_LITERAL(48, 492, 19), // "setZNearCoefficient"
QT_MOC_LITERAL(49, 512, 4), // "coef"
QT_MOC_LITERAL(50, 517, 23), // "setZClippingCoefficient"
QT_MOC_LITERAL(51, 541, 13), // "setOrthoZNear"
QT_MOC_LITERAL(52, 555, 1), // "z"
QT_MOC_LITERAL(53, 557, 10), // "orthoZNear"
QT_MOC_LITERAL(54, 568, 14), // "setSceneRadius"
QT_MOC_LITERAL(55, 583, 14), // "setSceneCenter"
QT_MOC_LITERAL(56, 598, 23), // "setSceneCenterFromPixel"
QT_MOC_LITERAL(57, 622, 19), // "setSceneBoundingBox"
QT_MOC_LITERAL(58, 642, 13), // "setPivotPoint"
QT_MOC_LITERAL(59, 656, 5), // "point"
QT_MOC_LITERAL(60, 662, 22), // "setPivotPointFromPixel"
QT_MOC_LITERAL(61, 685, 8), // "setFrame"
QT_MOC_LITERAL(62, 694, 28), // "ManipulatedCameraFrame*const"
QT_MOC_LITERAL(63, 723, 3), // "mcf"
QT_MOC_LITERAL(64, 727, 23), // "setKeyFrameInterpolator"
QT_MOC_LITERAL(65, 751, 1), // "i"
QT_MOC_LITERAL(66, 753, 26), // "KeyFrameInterpolator*const"
QT_MOC_LITERAL(67, 780, 3), // "kfi"
QT_MOC_LITERAL(68, 784, 17), // "addKeyFrameToPath"
QT_MOC_LITERAL(69, 802, 8), // "playPath"
QT_MOC_LITERAL(70, 811, 10), // "deletePath"
QT_MOC_LITERAL(71, 822, 9), // "resetPath"
QT_MOC_LITERAL(72, 832, 11), // "setFlySpeed"
QT_MOC_LITERAL(73, 844, 5), // "speed"
QT_MOC_LITERAL(74, 850, 15) // "onFrameModified"

    },
    "CGAL::qglviewer::Camera\0setPosition\0"
    "\0Vec\0pos\0setOrientation\0Quaternion\0q\0"
    "theta\0phi\0setUpVector\0up\0noMove\0"
    "setViewDirection\0direction\0lookAt\0"
    "target\0showEntireScene\0fitSphere\0"
    "center\0radius\0fitBoundingBox\0min\0max\0"
    "fitScreenRegion\0rectangle\0centerScene\0"
    "interpolateToZoomOnPixel\0pixel\0"
    "interpolateToFitScene\0interpolateTo\0"
    "Frame\0fr\0duration\0setType\0Type\0type\0"
    "setFieldOfView\0fov\0setHorizontalFieldOfView\0"
    "hfov\0setFOVToFitScene\0setAspectRatio\0"
    "aspect\0setScreenWidthAndHeight\0width\0"
    "height\0devicePixelRatio\0setZNearCoefficient\0"
    "coef\0setZClippingCoefficient\0setOrthoZNear\0"
    "z\0orthoZNear\0setSceneRadius\0setSceneCenter\0"
    "setSceneCenterFromPixel\0setSceneBoundingBox\0"
    "setPivotPoint\0point\0setPivotPointFromPixel\0"
    "setFrame\0ManipulatedCameraFrame*const\0"
    "mcf\0setKeyFrameInterpolator\0i\0"
    "KeyFrameInterpolator*const\0kfi\0"
    "addKeyFrameToPath\0playPath\0deletePath\0"
    "resetPath\0setFlySpeed\0speed\0onFrameModified"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_CGAL__qglviewer__Camera[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      40,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,  214,    2, 0x0a /* Public */,
       5,    1,  217,    2, 0x0a /* Public */,
       5,    2,  220,    2, 0x0a /* Public */,
      10,    2,  225,    2, 0x0a /* Public */,
      10,    1,  230,    2, 0x2a /* Public | MethodCloned */,
      13,    1,  233,    2, 0x0a /* Public */,
      15,    1,  236,    2, 0x0a /* Public */,
      17,    0,  239,    2, 0x0a /* Public */,
      18,    2,  240,    2, 0x0a /* Public */,
      21,    2,  245,    2, 0x0a /* Public */,
      24,    1,  250,    2, 0x0a /* Public */,
      26,    0,  253,    2, 0x0a /* Public */,
      27,    1,  254,    2, 0x0a /* Public */,
      29,    0,  257,    2, 0x0a /* Public */,
      30,    2,  258,    2, 0x0a /* Public */,
      34,    1,  263,    2, 0x0a /* Public */,
      37,    1,  266,    2, 0x0a /* Public */,
      39,    1,  269,    2, 0x0a /* Public */,
      41,    0,  272,    2, 0x0a /* Public */,
      42,    1,  273,    2, 0x0a /* Public */,
      44,    3,  276,    2, 0x0a /* Public */,
      44,    2,  283,    2, 0x2a /* Public | MethodCloned */,
      48,    1,  288,    2, 0x0a /* Public */,
      50,    1,  291,    2, 0x0a /* Public */,
      51,    1,  294,    2, 0x0a /* Public */,
      53,    0,  297,    2, 0x0a /* Public */,
      54,    1,  298,    2, 0x0a /* Public */,
      55,    1,  301,    2, 0x0a /* Public */,
      56,    1,  304,    2, 0x0a /* Public */,
      57,    2,  307,    2, 0x0a /* Public */,
      58,    1,  312,    2, 0x0a /* Public */,
      60,    1,  315,    2, 0x0a /* Public */,
      61,    1,  318,    2, 0x0a /* Public */,
      64,    2,  321,    2, 0x0a /* Public */,
      68,    1,  326,    2, 0x0a /* Public */,
      69,    1,  329,    2, 0x0a /* Public */,
      70,    1,  332,    2, 0x0a /* Public */,
      71,    1,  335,    2, 0x0a /* Public */,
      72,    1,  338,    2, 0x0a /* Public */,
      74,    0,  341,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, QMetaType::QReal, QMetaType::QReal,    8,    9,
    QMetaType::Void, 0x80000000 | 3, QMetaType::Bool,   11,   12,
    QMetaType::Void, 0x80000000 | 3,   11,
    QMetaType::Void, 0x80000000 | 3,   14,
    QMetaType::Void, 0x80000000 | 3,   16,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 3, QMetaType::QReal,   19,   20,
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3,   22,   23,
    QMetaType::Void, QMetaType::QRect,   25,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QPoint,   28,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 31, QMetaType::QReal,   32,   33,
    QMetaType::Void, 0x80000000 | 35,   36,
    QMetaType::Void, QMetaType::QReal,   38,
    QMetaType::Void, QMetaType::QReal,   40,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QReal,   43,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, QMetaType::QReal,   45,   46,   47,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   45,   46,
    QMetaType::Void, QMetaType::QReal,   49,
    QMetaType::Void, QMetaType::QReal,   49,
    QMetaType::Void, QMetaType::QReal,   52,
    QMetaType::QReal,
    QMetaType::Void, QMetaType::QReal,   20,
    QMetaType::Void, 0x80000000 | 3,   19,
    QMetaType::Bool, QMetaType::QPoint,   28,
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3,   22,   23,
    QMetaType::Void, 0x80000000 | 3,   59,
    QMetaType::Bool, QMetaType::QPoint,   28,
    QMetaType::Void, 0x80000000 | 62,   63,
    QMetaType::Void, QMetaType::UInt, 0x80000000 | 66,   65,   67,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::QReal,   73,
    QMetaType::Void,

       0        // eod
};

void CGAL::qglviewer::Camera::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Camera *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->setPosition((*reinterpret_cast< const Vec(*)>(_a[1]))); break;
        case 1: _t->setOrientation((*reinterpret_cast< const Quaternion(*)>(_a[1]))); break;
        case 2: _t->setOrientation((*reinterpret_cast< qreal(*)>(_a[1])),(*reinterpret_cast< qreal(*)>(_a[2]))); break;
        case 3: _t->setUpVector((*reinterpret_cast< const Vec(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 4: _t->setUpVector((*reinterpret_cast< const Vec(*)>(_a[1]))); break;
        case 5: _t->setViewDirection((*reinterpret_cast< const Vec(*)>(_a[1]))); break;
        case 6: _t->lookAt((*reinterpret_cast< const Vec(*)>(_a[1]))); break;
        case 7: _t->showEntireScene(); break;
        case 8: _t->fitSphere((*reinterpret_cast< const Vec(*)>(_a[1])),(*reinterpret_cast< qreal(*)>(_a[2]))); break;
        case 9: _t->fitBoundingBox((*reinterpret_cast< const Vec(*)>(_a[1])),(*reinterpret_cast< const Vec(*)>(_a[2]))); break;
        case 10: _t->fitScreenRegion((*reinterpret_cast< const QRect(*)>(_a[1]))); break;
        case 11: _t->centerScene(); break;
        case 12: _t->interpolateToZoomOnPixel((*reinterpret_cast< const QPoint(*)>(_a[1]))); break;
        case 13: _t->interpolateToFitScene(); break;
        case 14: _t->interpolateTo((*reinterpret_cast< const Frame(*)>(_a[1])),(*reinterpret_cast< qreal(*)>(_a[2]))); break;
        case 15: _t->setType((*reinterpret_cast< Type(*)>(_a[1]))); break;
        case 16: _t->setFieldOfView((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 17: _t->setHorizontalFieldOfView((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 18: _t->setFOVToFitScene(); break;
        case 19: _t->setAspectRatio((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 20: _t->setScreenWidthAndHeight((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< qreal(*)>(_a[3]))); break;
        case 21: _t->setScreenWidthAndHeight((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 22: _t->setZNearCoefficient((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 23: _t->setZClippingCoefficient((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 24: _t->setOrthoZNear((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 25: { qreal _r = _t->orthoZNear();
            if (_a[0]) *reinterpret_cast< qreal*>(_a[0]) = std::move(_r); }  break;
        case 26: _t->setSceneRadius((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 27: _t->setSceneCenter((*reinterpret_cast< const Vec(*)>(_a[1]))); break;
        case 28: { bool _r = _t->setSceneCenterFromPixel((*reinterpret_cast< const QPoint(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        case 29: _t->setSceneBoundingBox((*reinterpret_cast< const Vec(*)>(_a[1])),(*reinterpret_cast< const Vec(*)>(_a[2]))); break;
        case 30: _t->setPivotPoint((*reinterpret_cast< const Vec(*)>(_a[1]))); break;
        case 31: { bool _r = _t->setPivotPointFromPixel((*reinterpret_cast< const QPoint(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        case 32: _t->setFrame((*reinterpret_cast< ManipulatedCameraFrame*const(*)>(_a[1]))); break;
        case 33: _t->setKeyFrameInterpolator((*reinterpret_cast< uint(*)>(_a[1])),(*reinterpret_cast< KeyFrameInterpolator*const(*)>(_a[2]))); break;
        case 34: _t->addKeyFrameToPath((*reinterpret_cast< uint(*)>(_a[1]))); break;
        case 35: _t->playPath((*reinterpret_cast< uint(*)>(_a[1]))); break;
        case 36: _t->deletePath((*reinterpret_cast< uint(*)>(_a[1]))); break;
        case 37: _t->resetPath((*reinterpret_cast< uint(*)>(_a[1]))); break;
        case 38: _t->setFlySpeed((*reinterpret_cast< qreal(*)>(_a[1]))); break;
        case 39: _t->onFrameModified(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject CGAL::qglviewer::Camera::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_CGAL__qglviewer__Camera.data,
    qt_meta_data_CGAL__qglviewer__Camera,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *CGAL::qglviewer::Camera::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *CGAL::qglviewer::Camera::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_CGAL__qglviewer__Camera.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int CGAL::qglviewer::Camera::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 40)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 40;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 40)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 40;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
