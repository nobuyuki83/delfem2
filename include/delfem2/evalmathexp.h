/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


#ifndef DFM2_EVALMATHEXP_H
#define DFM2_EVALMATHEXP_H

#include <iostream>
#include <string>
#include <vector>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

namespace evalmathexp {

class CCmd {
public:
  virtual ~CCmd() = default;
  virtual bool DoOperation(std::vector<double> &stack) = 0;
  virtual void SetValue(const double &val) = 0;
};

class CKey
{	// class for storing values of key
public:
  CKey(const std::string& name, double val)
      : m_Name(name), m_Val(val){}
  std::string m_Name;
  std::vector<unsigned int> m_aiCmd;
  double m_Val;
};

}

// ----------------------

class CMathExpressionEvaluator
{
public:
	CMathExpressionEvaluator(){ m_is_valid = false; }
	CMathExpressionEvaluator(const std::string& exp){ this->SetExp(exp); }
	~CMathExpressionEvaluator();

	// --------------------
	bool SetExp(const std::string& key_name);
	void SetKey(const std::string& key_name, double val);
	bool IsKeyUsed(const std::string& key_name){
		for(unsigned int ikey=0;ikey<m_aKey.size();ikey++){
			if( m_aKey[ikey].m_Name == key_name ){
				if( !m_aKey[ikey].m_aiCmd.empty() ) return true;
				return false; 
			}
		}
		return false;
	}
	double Eval() const;
private:
	bool m_is_valid;
	std::string m_sExp;
	std::vector<evalmathexp::CCmd*> m_apCmd;	// 逆ポーランド記法された演算子や値の列
	std::vector<evalmathexp::CKey> m_aKey;	// 文字列の名前と値が，どのIndexのコマンドに格納されているか
};
  
}

#ifndef DFM2_STATIC_LIBRARY
#  include "delfem2/evalmathexp.cpp"
#endif

#endif	// !defind EVALMATHEXP
