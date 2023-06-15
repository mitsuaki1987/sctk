.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

インストール方法
================
 
要件
----

使用要件はオリジナルの QuantumESPRESSO_ と同じである。
          
インストール手順
----------------

#. 以下のようにして、SCTKのリポジトリからソースコードをgitクローンする。
   これは QuantumESPRESSO_ 本体からフォークされたものである。

   .. code-block:: bash

      $ git clone git://git.osdn.net/gitroot/sctk/sctk.git

#. オリジナルの QuantumESPRESSO_ の場合と同様に ``configure`` で環境設定をする。
               
   .. code-block:: bash

       $ ./configure
       
#. メイク

   .. code-block:: bash

       $ make sctk
               
   実行ファイル名は ``sctk.x`` である。

