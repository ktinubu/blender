/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __BVM_CODEGEN_H__
#define __BVM_CODEGEN_H__

/** \file bvm_codegen.h
 *  \ingroup bvm
 */

#include <set>
#include <vector>

#include "MEM_guardedalloc.h"

#include "compiler.h"
#include "node_graph.h"

#include "bvm_instruction_list.h"

#include "util_opcode.h"
#include "util_string.h"

namespace blenvm {

struct FunctionBVM;
struct NodeGraph;
struct NodeInstance;
struct TypeDesc;

struct Value {
	Value(StackIndex m_stack_index);
	~Value();
	
	StackIndex stack_index() const { return m_stack_index; }
	
	const NodeConstant *constant_value() const { return m_constant_value; }
	void set_constant_value(const NodeConstant *value);
	
private:
	StackIndex m_stack_index;
	NodeConstant *m_constant_value;
};

struct BVMCodeGenerator : public CodeGenerator {
	typedef std::vector<Value> ValueStack;
	typedef std::map<ValueHandle, Value*> HandleValueMap;
	typedef std::vector<int> StackUsers;
	typedef std::vector<Value*> Arguments;
	
	BVMCodeGenerator();
	~BVMCodeGenerator();
	
	void finalize_function();
	void debug_function(FILE *file);
	
	void node_graph_begin(const string &name, const NodeGraph *graph, bool use_globals);
	void node_graph_end();
	
	void store_return_value(size_t output_index, const TypeSpec *typespec, ValueHandle value);
	ValueHandle map_argument(size_t input_index, const TypeSpec *typespec);
	
	ValueHandle alloc_node_value(const TypeSpec *typespec, const string &name);
	ValueHandle create_constant(const TypeSpec *typespec, const NodeConstant *node_value);
	
	void eval_node(const NodeType *nodetype,
	               ArrayRef<ValueHandle> input_args,
	               ArrayRef<ValueHandle> output_args);
	
protected:
	static ValueHandle get_handle(const Value *value);
	
	StackIndex find_stack_index(int size) const;
	StackIndex assign_stack_index(const TypeDesc &typedesc);
	Value *create_value(const TypeSpec *typespec);
	
private:
	ValueStack m_values;
	HandleValueMap m_value_map;
	StackUsers m_stack_users;
	Arguments m_input_args;
	Arguments m_output_args;
	
	InstructionList m_instructions;
};


/* ========================================================================= */


typedef std::map<ConstInputKey, StackIndex> InputIndexMap;
typedef std::map<ConstOutputKey, StackIndex> OutputIndexMap;
typedef std::map<ConstOutputKey, int> OutputUsersMap;
typedef std::map<const NodeInstance *, OutputSet> NodeDependencyMap;
typedef std::map<const NodeBlock *, OutputSet> BlockDependencyMap;

typedef std::vector<int> StackUsers;

struct BVMCompilerBase {
	struct BlockInfo {
		int entry_point;
	};
	typedef std::map<const NodeBlock*, BlockInfo> BlockInfoMap;
	
	BVMCompilerBase();
	virtual ~BVMCompilerBase();
	
protected:
	void calc_node_dependencies(const NodeInstance *node, BlockDependencyMap &block_deps_map);
	void calc_block_dependencies(const NodeBlock *block, BlockDependencyMap &block_deps_map);
	void calc_symbol_scope(const NodeGraph &graph);
	
	StackIndex find_stack_index(int size) const;
	StackIndex assign_stack_index(const TypeDesc &typedesc);
	void get_local_arg_indices(const NodeInstance *node, const NodeBlock *local_block);
	void resolve_node_block_symbols(const NodeBlock *block);
	void resolve_symbols(const NodeGraph &graph);
	
	virtual void push_opcode(OpCode op) const = 0;
	virtual void push_stack_index(StackIndex arg) const = 0;
	virtual void push_jump_address(int address) const = 0;
	
	virtual void push_float(float f) const = 0;
	virtual void push_float3(float3 f) const = 0;
	virtual void push_float4(float4 f) const = 0;
	virtual void push_int(int i) const = 0;
	virtual void push_matrix44(matrix44 m) const = 0;
	virtual void push_string(const char *s) const = 0;
	
	virtual int current_address() const = 0;
	
	void push_constant(const NodeConstant *value) const;
	
	void codegen_value(const NodeConstant *value, StackIndex offset) const;
	int codegen_node_block(const NodeBlock &block);
	int codegen_graph(const NodeGraph &graph);
	
protected:
	BlockInfoMap block_info;
	NodeDependencyMap node_deps_map;
	OutputUsersMap output_users;
	
	InputIndexMap input_index;
	OutputIndexMap output_index;
	StackUsers stack_users;

	MEM_CXX_CLASS_ALLOC_FUNCS("BVM:BVMCompiler")
};

struct BVMCompiler : public BVMCompilerBase {
	BVMCompiler();
	~BVMCompiler();
	
	FunctionBVM *compile_function(const NodeGraph &graph);
	
protected:
	void push_opcode(OpCode op) const;
	void push_stack_index(StackIndex arg) const;
	void push_jump_address(int address) const;
	
	void push_float(float f) const;
	void push_float3(float3 f) const;
	void push_float4(float4 f) const;
	void push_int(int i) const;
	void push_matrix44(matrix44 m) const;
	void push_string(const char *s) const;
	
	int current_address() const;
	
private:
	FunctionBVM *fn;
};

struct DebugGraphvizBVMCompiler : public BVMCompilerBase {
	DebugGraphvizBVMCompiler();
	~DebugGraphvizBVMCompiler();
	
	void compile_function(const NodeGraph &graph, FILE *file, const string &label);
	
protected:
	void push_opcode(OpCode op) const;
	void push_stack_index(StackIndex arg) const;
	void push_jump_address(int address) const;
	
	void push_float(float m_file) const;
	void push_float3(float3 m_file) const;
	void push_float4(float4 m_file) const;
	void push_int(int i) const;
	void push_matrix44(matrix44 m) const;
	void push_string(const char *s) const;
	
	int current_address() const;
	
	const char *get_arg_name() const;
	bool is_arg_output() const;
	
	void init_graph(const string &label);
	void close_graph();
	
	void init_node();
	void close_node();
	
private:
	FILE *m_file;
	mutable int m_current_address;
	mutable const NodeType *m_current_opnode;
	mutable int m_current_arg;
};

} /* namespace blenvm */

#endif /* __BVM_CODEGEN_H__ */