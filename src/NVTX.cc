#include "NVTX.H"

#if HAVE_LIBCUDA==1

const uint32_t colors[] = { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000, 0x00ffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

unsigned nvTracer::count = 0;

nvTracer::nvTracer(const char *name)
{
  int color_id = count++ % num_colors;

  nvtxEventAttributes_t eventAttrib = {0};

  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = colors[color_id];
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii = name;

  nvtxRangePushEx(&eventAttrib);
}

#endif
