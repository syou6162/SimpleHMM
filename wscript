VERSION = '0.0.1'
APPNAME = 'simple_hmm'

srcdir = '.'
blddir = 'build'

def set_options(ctx):
  ctx.tool_options('compiler_cxx')

def configure(ctx):
  ctx.check_tool('compiler_cxx')
  ctx.env.append_value('std', ['c++0x'])
  ctx.check_cxx(lib = 'glog')
  ctx.check_cxx(lib = 'gflags')
  ctx.env.CXXFLAGS += ['-O3', '-Wall', '-W']

def build(bld):
  bld.recurse('src')
