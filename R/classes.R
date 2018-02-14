require(methods)

# S4 classes

setClass("DataList",
representation("list")
)

# Set inheritance
# The LargeDataObject class is set in limma and provides a show method

setIs("DataList", "LargeDataObject")
